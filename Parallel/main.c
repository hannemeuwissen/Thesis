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

        if(floor(M/nprocs) < degree){
            printf("Invalid input: the dimensions must define a tall skinny matrix on every process (dimension on process: %d x %d).\n", M/nprocs, degree);
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

    // /* Generate part of transition matrix for calling process */
    // sparse_CSR A = generate_regular_graph_part_csr(m, M, nnz);

    /* Test: read from file (each process)*/
    sparse_CSR A;
    read_CSR(&A, "smallcsr.txt");

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
    
    /* Normalize start vector */
    double local_dot = cblas_ddot(m, v, 1, v, 1);
    double global_dot;
    MPI_Allreduce(&local_dot, &global_dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_norm = sqrt(global_dot);
    for(int i=0;i<m;i++){v[i] /= global_norm;}

    /* CA-Arnoldi(s, steps) (note: no restarting, final degree = s*steps) */
    for(int block = 0;block < steps;block++){
        printf("** Block %d\n", block);

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

            // if(!myid){
            //     print_matrix(mathcalQ, (steps*s + 1), m);
            // }

            /* Save last vector in v */
            memcpy(v, V + s*m, m*sizeof(double));
            free(V);

            /* Calculate mathcal H (note: only process 0 calculates H, and final H is not transposed!)*/
            if(!myid){

                // print_matrix(R_, s+1, s+1);
                mathcalH = malloc((s+1)*s*sizeof(double));
                double * R = malloc(s*s*sizeof(double)); 
                get_R(R, R_, s+1);
                double * B_ = malloc(s*(s+1)*sizeof(double));
                set_B_(B_, s);
                calc_hess_on_transpose(mathcalH, R_, B_, R, s+1, s);
                free(R);
                free(B_);

                // print_matrix(mathcalH, s+1, s);
            }
            free(R_);
        }else{
            /* Matrix powers kernel (note: saved as transpose - vectors in rows!) */
            V = malloc(s*m*sizeof(double));
            matrix_powers(A, v, V, s, m, myid, nprocs, start, end, MPI_COMM_WORLD);

            /* Block-CGS to orthogonalize compared to previous blocks */
            mathcalR_ = malloc((block*s + 1)*s*sizeof(double));
            bgs_on_transpose(mathcalQ, V, mathcalR_, m, block*s + 1, s, MPI_COMM_WORLD); // note: mathcal R is not transposed
            // MPI_Barrier(MPI_COMM_WORLD);
            printf("Process %d finished bgs\n", myid);
            
            /* Orthogonalize block using parallel CA-TSQR */
            R_ = malloc(s*s*sizeof(double)); 
            TSQR_on_transpose(V, m, s, R_, myid, nprocs, MPI_COMM_WORLD); // note: resulting R is transposed!
            // print_matrix(R_, s, s);
            // MPI_Barrier(MPI_COMM_WORLD); 

            /* Set mathcal Q (note: saved as transpose - vectors in rows!) */
            memcpy(mathcalQ + (1 + block*s)*m, V, s*m*sizeof(double));

            /* Save last vector in v */
            memcpy(v, V + (s-1)*m, m*sizeof(double));
            free(V);

            /* Update mathcal H */
            if(!myid){
                // print_matrix(R_, s, s);
                update_hess_on_transpose(&mathcalH, mathcalR_, R_, s, block);
            }
            free(R_);
            free(mathcalR_);
        }

        // printf("Process %d finished block %d\n", myid, block);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    if(!myid){
        print_matrix(mathcalQ, (steps*s + 1), m);
        print_matrix(mathcalH, (steps*s + 1), steps*s);
    }

    free(mathcalQ);
    free(mathcalH);

    free(start);
    free(end);
    
    MPI_Finalize();
    return 0;
}