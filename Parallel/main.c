/**
 * @file main.c
 * @author Hanne Meuwissen (22307813)
 * @brief Code for thesis at Trinity College Dublin.
 * @version 1.0
 * @date 2023-06-02
 */
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
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
    int degree,s,M,nnz;
    char filename_v[100];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(!myid){ /* Read input */
        float logprocs = log2(nprocs);
        if(ceil(logprocs) != floor(logprocs)){
            fprintf(stderr,"Error: The number of processes needs to be a power of 2 (because of TSQR).\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        parse_command_line_regular(argc, argv, &M, &nnz, filename_v, &degree, &s, MPI_COMM_WORLD);
        if(degree%s != 0){
            printf("Invalid input: the degree of the Krylov subspace should be a multiple of the blocksize (s)\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if((M<=0) || (degree<=0) || (M<degree)){
            printf("Invalid input: the dimensions must define a tall skinny matrix.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if(floor(M/nprocs) < s+1){
            printf("Invalid input: the dimensions must define a tall skinny matrix on every process (dimension on process in step 0: %d x %d).\nSuggestion: lower the blocksize.\n", M/nprocs, s+1);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    /* Broadcast input to all processes */
    MPI_Bcast(&degree, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(filename_v, 100, MPI_CHAR, 0, MPI_COMM_WORLD);

    /* Determine the start index and size of part for calling process */
    int * start = malloc(nprocs*sizeof(int));
    int * end = malloc(nprocs*sizeof(int));
    get_indices(M, nprocs, start, end);
    int m = end[myid] - start[myid] + 1;

    /* Generate part of transition matrix for calling process */
    sparse_CSR A = generate_regular_graph_part_csr(m, M, nnz);

    // /* Test: read from file (each process)*/
    // sparse_CSR A;
    // read_CSR(&A, "smallcsr4.txt");

    // if(!myid){
    //     print_CSR(&A);
    // }

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
    
    /* Normalize start vector */
    double local_dot = cblas_ddot(m, v, 1, v, 1);
    double global_dot;
    MPI_Allreduce(&local_dot, &global_dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_norm = sqrt(global_dot);
    for(int i=0;i<m;i++){v[i] /= global_norm;}

    /* CA-Arnoldi(s, steps) (note: no restarting, final degree = s*steps) */
    for(int block = 0;block < steps;block++){
        // printf("** Block %d\n", block);

        if(!block){
            /* Matrix powers kernel (note: saved as transpose - vectors in rows!)*/
            V = malloc((s+1)*m*sizeof(double));
            memcpy(V, v, m*sizeof(double));
            matrix_powers(A, v, V + m, s, m, myid, nprocs, start, end, MPI_COMM_WORLD);

            /* Orthogonalize first block using parallel CA-TSQR */
            R_ = malloc((s+1)*(s+1)*sizeof(double)); 
            TSQR_on_transpose(V, m, s + 1, R_, myid, nprocs, MPI_COMM_WORLD); // note: resulting R is transposed!
            // MPI_Barrier(MPI_COMM_WORLD);

            /* Set mathcal Q (note: saved as transpose - vectors in rows!) */
            memcpy(mathcalQ, V, (s+1)*m*sizeof(double));

            /* Save last vector in v */
            memcpy(v, V + s*m, m*sizeof(double));
            free(V);

            /* Calculate mathcal H (note: only process 0 calculates H, and final H is not transposed!)*/
            if(!myid){
                mathcalH = malloc((s+1)*s*sizeof(double));
                double * R = malloc(s*s*sizeof(double)); 
                get_R(R, R_, s+1);
                double * B_ = malloc(s*(s+1)*sizeof(double));
                set_B_(B_, s);
                calc_hess_on_transpose(mathcalH, R_, B_, R, s+1, s);
                free(R);
                free(B_);

                print_matrix(mathcalH, s+1, s);

                breakdown = breakdown_check(mathcalH, s, block, tol);
                // printf("Breakdown block %d: %d\n", block, breakdown);
            }
            MPI_Bcast(&breakdown, 1, MPI_INT, 0, MPI_COMM_WORLD);
            free(R_);
        }else{
            /* Matrix powers kernel (note: saved as transpose - vectors in rows!) */
            V = malloc(s*m*sizeof(double));
            matrix_powers(A, v, V, s, m, myid, nprocs, start, end, MPI_COMM_WORLD);

            /* Block-CGS to orthogonalize compared to previous blocks */
            mathcalR_ = malloc((block*s + 1)*s*sizeof(double));
            bgs_on_transpose(mathcalQ, V, mathcalR_, m, block*s + 1, s, MPI_COMM_WORLD); // note: mathcal R is not transposed
            // MPI_Barrier(MPI_COMM_WORLD);
            
            /* Orthogonalize block using parallel CA-TSQR */
            R_ = malloc(s*s*sizeof(double)); 
            TSQR_on_transpose(V, m, s, R_, myid, nprocs, MPI_COMM_WORLD); // note: resulting R is transposed!
            // MPI_Barrier(MPI_COMM_WORLD); 

            /* Set mathcal Q (note: saved as transpose - vectors in rows!) */
            memcpy(mathcalQ + (1 + block*s)*m, V, s*m*sizeof(double));

            /* Save last vector in v */
            memcpy(v, V + (s-1)*m, m*sizeof(double));
            free(V);

            /* Update mathcal H */
            if(!myid){
                update_hess_on_transpose(&mathcalH, mathcalR_, R_, s, block);
                print_matrix(mathcalH, s*(block + 1) + 1, s*(block + 1));
                breakdown = breakdown_check(mathcalH, s, block, tol);
                // printf("Breakdown block %d: %d\n", block, breakdown);
            }
            MPI_Bcast(&breakdown, 1, MPI_INT, 0, MPI_COMM_WORLD);

            free(R_);
            free(mathcalR_);

        }

        MPI_Barrier(MPI_COMM_WORLD);

        if(breakdown != -1){
            /* Set matching vectors in mathcalQ to zero */
            printf("set from column %d (%d elements) to zero\n", breakdown, ((1 + degree) - (breakdown -1))*m);
            // memset(mathcalQ + breakdown*m, 0, ((1 + degree) - (breakdown -1))*m*sizeof(double));

            // if(!myid){
            //     /* Fill mathcall H with zeros */
            //     double * temp = malloc((degree + 1)*degree*sizeof(double));
            //     for(int i=0;i<=breakdown;i++){
            //         memcpy(temp + i*degree, mathcalH + i*(s*(block+1)), (s*(block+1))*sizeof(double));
            //         memset(temp + (s*(block+1)) + i*degree, 0, (degree - (breakdown - 1))*sizeof(double));
            //     }
            //     memset(temp + (breakdown+1)*degree, 0, degree*((1 + degree) - (breakdown -1))*sizeof(double));
            //     free(mathcalH);
            //     mathcalH = temp;
            // }
            printf("Breakdown row %d in block %d\n", breakdown, block);

            break;
        }
    }

    if(!myid){ /* Print out results */
        printf("Part of Q process 0:\n");
        print_matrix_transposed(mathcalQ, (steps*s + 1), m);
        printf("Hessenberg:\n");
        print_matrix(mathcalH, (steps*s + 1), steps*s);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // if(myid == 1){ /* Print out results */
    //     printf("Part of Q process 1:\n");
    //     print_matrix_transposed(mathcalQ, (steps*s + 1), m);
    // }


    free(mathcalQ);
    if(!myid){free(mathcalH);}
    free(start);
    free(end);
    
    MPI_Finalize();
    return 0;
}