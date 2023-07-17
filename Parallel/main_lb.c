/**
 * @file main_lb.c
 * @brief Main code to investigate CA-Arnoldi with load-balancing, part of Thesis project in 
 * High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 1.0
 * @date 2023-07-13
 */
#include<stdlib.h>
#include<stdio.h>
#include<mpi.h>
#include<string.h>
#include<math.h>
#include<mkl_cblas.h>
#include<mkl_types.h>
#include"parse.h"
#include"matrix.h"
#include"graph.h"
#include"sparse.h"
#include"tsqr_mpi.h"
#include"bgs.h"
#include"hess.h"
#include"mkl.h"

int main(int argc, char **argv){
    const double tol = 1.0E-10;  
    int myid, nprocs;
    int degree,original_degree,s,M,total_nnz,lb;
    char filename_v[100], filename_A[100];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(!myid){ /* Read input */
        float logprocs = log2(nprocs);
        if(ceil(logprocs) != floor(logprocs)){
            fprintf(stderr,"Error: The number of processes needs to be a power of 2 (because of TSQR).\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        parse_command_line_lb(argc, argv, filename_A, filename_v, &degree, &s, &lb, MPI_COMM_WORLD);
        if((degree < 1) || (s > degree)){
            printf("Invalid input: the degree of the Krylov subspace should be at least 1 and the blocksize should be smaller\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        /* Handle the case if s doesn't divide degree */
        original_degree = degree;
        if(degree%s != 0){
            while(degree%s != 0){
                degree++;
            }
        }
    }

    /* Broadcast input to all processes */
    MPI_Bcast(&degree, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&original_degree, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(filename_v, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(filename_A, 100, MPI_CHAR, 0, MPI_COMM_WORLD);

    /* Read first CSR data from the file */
    sparse_CSR A;
    read_CSR_data(&A, filename_A);
    M = A.nrows;
    if(!myid){
        if(M != A.ncols){
            printf("Invalid input file: The matrix doesn't represent a graph (non-square).\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if((M<=0) || (M<original_degree)){
            printf("Invalid input: the dimensions must define a tall skinny matrix.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    /* Determine the start index, end index and size of part for calling process */
    int * start = malloc(nprocs*sizeof(int));
    int * end = malloc(nprocs*sizeof(int));
    if(!lb){
        get_indices(A.nrows, nprocs, start, end);
    }else{
        get_indices_load_balanced(A, nprocs, start, end);
    }
    int m = end[myid] - start[myid] + 1;
    if(m < s+1){
        printf("Invalid input: the dimensions must define a tall skinny matrix on every process (dimension on process in step 0: %d x %d).\nSuggestion: lower the blocksize.\n", M/nprocs, s+1);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* Read part of transition matrix for calling process */
    sparse_CSR A;
    read_CSR_part(&A, filename_A, start[myid], end[myid]);
    printf("Process %d is done reading its part!\n", myid);

    /* Initialize arrays */
    int steps = degree/s;
    double *V;
    double *R_;
    double *mathcalQ = malloc((1 + steps*s)*m*sizeof(double));
    double *mathcalH;
    double *mathcalR_;
    double *v = malloc(m*sizeof(double));
    read_matrix_from_file(filename_v, start[myid], v, m, 1);

    /* Initialize breakdown indicator */
    int breakdown = -1;

    /* Initialize timing variables */
    double * mp_times = malloc(steps*sizeof(double));
    double * bgs_times = malloc((steps-1)*sizeof(double));
    double * tsqr_times = malloc(steps*sizeof(double));
    double * hess_times;
    if(!myid){hess_times = malloc(steps*sizeof(double));}
    double t1 = MPI_Wtime();
    double tbeg, tend;
    
    /* Normalize start vector */
    double local_dot = cblas_ddot(m, v, 1, v, 1);
    double global_dot;
    MPI_Allreduce(&local_dot, &global_dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_norm = sqrt(global_dot);
    for(int i=0;i<m;i++){v[i] /= global_norm;}

    /* CA-Arnoldi(s, steps) (note: no restarting, final degree = s*steps) */
    int block = 0;
    for(block = 0;block < steps;block++){
        // printf("** Block %d\n", block);

        if(!block){

            /* Matrix powers kernel (note: saved as transpose - vectors in rows!)*/
            tbeg = MPI_Wtime();
            V = malloc((s+1)*m*sizeof(double));
            memcpy(V, v, m*sizeof(double));
            matrix_powers(A, v, V + m, s, m, myid, nprocs, start, end, MPI_COMM_WORLD);
            tend = MPI_Wtime();
            mp_times[block] = tend - tbeg;

            /* Orthogonalize first block using parallel CA-TSQR */
            tbeg = MPI_Wtime();
            R_ = malloc((s+1)*(s+1)*sizeof(double)); 
            TSQR_on_transpose(V, m, s + 1, R_, myid, nprocs, MPI_COMM_WORLD); // note: resulting R is transposed!
            tend = MPI_Wtime();
            tsqr_times[block] = tend - tbeg;

            /* Set mathcal Q (note: saved as transpose - vectors in rows!) */
            memcpy(mathcalQ, V, (s+1)*m*sizeof(double));

            /* Save last vector in v */
            memcpy(v, V + s*m, m*sizeof(double));
            free(V);

            /* Calculate mathcal H (note: only process 0 calculates H, and final H is not transposed!)*/
            if(!myid){
                tbeg = MPI_Wtime();
                mathcalH = malloc((s+1)*s*sizeof(double));
                double * R = malloc(s*s*sizeof(double)); 
                get_R(R, R_, s+1);
                double * B_ = malloc(s*(s+1)*sizeof(double));
                set_B_(B_, s);
                calc_hess_on_transpose(mathcalH, R_, B_, R, s+1, s);
                free(R);
                free(B_);
                breakdown = breakdown_check(mathcalH, s, block, tol);
                tend = MPI_Wtime();
                hess_times[block] = tend-tbeg;
            }

            MPI_Bcast(&breakdown, 1, MPI_INT, 0, MPI_COMM_WORLD);
            free(R_);

        }else{
            /* Matrix powers kernel (note: saved as transpose - vectors in rows!) */
            tbeg = MPI_Wtime();
            V = malloc(s*m*sizeof(double));
            matrix_powers(A, v, V, s, m, myid, nprocs, start, end, MPI_COMM_WORLD);
            tend = MPI_Wtime();
            mp_times[block] = tend-tbeg;

            /* Block-CGS to orthogonalize compared to previous blocks */
            tbeg = MPI_Wtime();
            mathcalR_ = malloc((block*s + 1)*s*sizeof(double));
            bgs_on_transpose(mathcalQ, V, mathcalR_, m, block*s + 1, s, MPI_COMM_WORLD); // note: mathcal R is not transposed
            tend = MPI_Wtime();
            bgs_times[block-1] = tend-tbeg;
            
            /* Orthogonalize block using parallel CA-TSQR */
            tbeg = MPI_Wtime();
            R_ = malloc(s*s*sizeof(double)); 
            TSQR_on_transpose(V, m, s, R_, myid, nprocs, MPI_COMM_WORLD); // note: resulting R is transposed!
            tend = MPI_Wtime(); 
            tsqr_times[block] = tend - tbeg;

            /* Set mathcal Q (note: saved as transpose - vectors in rows!) */
            memcpy(mathcalQ + (1 + block*s)*m, V, s*m*sizeof(double));

            /* Save last vector in v */
            memcpy(v, V + (s-1)*m, m*sizeof(double));
            free(V);

            /* Update mathcal H */
            if(!myid){
                tbeg = MPI_Wtime();
                update_hess_on_transpose(&mathcalH, mathcalR_, R_, s, block);
                breakdown = breakdown_check(mathcalH, s, block, tol);
                tend = MPI_Wtime();
                hess_times[block] = tend-tbeg;
            }
            MPI_Bcast(&breakdown, 1, MPI_INT, 0, MPI_COMM_WORLD);

            free(R_);
            free(mathcalR_);

        }

        MPI_Barrier(MPI_COMM_WORLD);

        /* Handle extra data in H and Q in case of breakdown */
        if(breakdown != -1){
            /* Set matching vectors in mathcalQ to zero */
            memset(mathcalQ + breakdown*m, 0, ((1 + degree) - breakdown)*m*sizeof(double));

            if(!myid){
                /* Fill mathcall H with zeros */
                double * temp = malloc((degree + 1)*degree*sizeof(double));
                for(int i=0;i<breakdown;i++){
                    memcpy(temp + i*degree, mathcalH + i*(s*(block+1)), (s*(block+1))*sizeof(double));
                    memset(temp + breakdown + i*degree, 0, (degree - breakdown)*sizeof(double));
                }
                memset(temp + breakdown*degree, 0, degree*((1 + degree) - breakdown)*sizeof(double));
                free(mathcalH);
                mathcalH = temp;
            }

            break;
        }
    }

    /* Remove extra data in H and Q in case s doesn't divide degree */
    if(original_degree != degree){
        /* Remove necessary columns from mathcalQ */
        double * temp = malloc((original_degree + 1)*m*sizeof(double));
        memcpy(temp, mathcalQ, (original_degree+1)*m*sizeof(double));
        free(mathcalQ);
        mathcalQ = temp;

        if(!myid){
            /* Remove necessary rows and columns from mathcalH */
            double * temp = malloc((original_degree + 1)*original_degree*sizeof(double));
            for(int i=0;i<(original_degree+1);i++){
                memcpy(temp + i*original_degree, mathcalH + i*degree, original_degree*sizeof(double));
            }
            free(mathcalH);
            mathcalH = temp;
        }
    }

    double t2 = MPI_Wtime();

    if(!myid){ /* Print out results */
        // printf("Part of Q process 0:\n");
        // print_matrix_transposed(mathcalQ, original_degree +1, m);
        printf("Total runtime process %d: %lf\n", myid, t2 - t1);
        printf("Average time matrix powers process %d: %lf\n", myid, average(mp_times, block));
        printf("Average time block (classical) Gram-Schmidt process %d: %lf\n", myid, average(bgs_times, block - 1));
        printf("Average time TSQR process %d: %lf\n", myid, average(tsqr_times, block));
        printf("Upper Hessenberg: %lf\n", average(hess_times, block - 1));
        printf("Hessenberg:\n");
        print_matrix(mathcalH, original_degree + 1, original_degree);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf("Total runtime process %d: %lf\n", myid, t2 - t1);
    printf("Average time matrix powers process %d: %lf\n", myid, average(mp_times, block));
    printf("Average time block (classical) Gram-Schmidt process %d: %lf\n", myid, average(bgs_times, block - 1));
    printf("Average time TSQR process %d: %lf\n", myid, average(tsqr_times, block));

    if(myid == 1){ /* Print out results */
        // printf("Part of Q process 1:\n");
        // print_matrix_transposed(mathcalQ, original_degree + 1, m);
        // printf("Runtime process %d: %lf\n", myid, t2 - t1);
    }

    free(mathcalQ);
    if(!myid){free(mathcalH);}
    free(start);
    free(end);
    
    MPI_Finalize();
    return 0;
}