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
    int degree,s,M,N,nnz;
    char filename_v[100];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(!myid){
        float logprocs = log2(nprocs);
        if(ceil(logprocs) != floor(logprocs)){
            fprintf(stderr,"Error: The number of processes needs to be a power of 2 (because of TSQR).\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        parse_command_line_regular(argc, argv, &M, &N, &nnz, filename_v, &degree, &s, MPI_COMM_WORLD);
        if(degree%s != 0){
            printf("Invalid input: the degree of the Krylov subspace should be a multiple of the blocksize (s)\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if((M<=0) || (N<=0) || (M<N)){
            printf("Invalid input: the dimensions must define a tall skinny matrix.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if(floor(M/nprocs) < N){
            printf("Invalid input: the dimensions must define a tall skinny matrix on every process (dimension on process: %d x %d).\n", M/nprocs, N);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Bcast(&degree, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(filename_v, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
    int steps = degree/s;
    int start, end;
    decomp1d(M, nprocs, myid, &start, &end);
    int m = end - start + 1;

    /* Generate part of transition matrix for calling process */
    sparse_CSR A = generate_regular_graph_part_csr(m, M, nnz);

    /* Initialize arrays */
    double *V;
    double *R_;
    double *mathcalQ = malloc((1 + steps*s)*m*sizeof(double));
    double *mathcalH;
    double *mathcalR_;
    double *v = malloc(m*sizeof(double));
    read_matrix_from_file(filename_v, start, v, m, 1);
    // NOTE: does this vector need to be normalized? -> don't thinks so because normalized by TSQR?

    /* CA-Arnoldi(s, steps) (note: no restarting, final degree = s*steps) */
    for(int block = 0;block < steps;block++){

        if(!block){
            /* Matrix powers kernel (note: saved as transpose - vectors in rows!)*/
            V = malloc((s+1)*m*sizeof(double));
            memcpy(V, v, m*sizeof(double));
            matrix_powers(A, v, V + m, s, m, myid, nprocs, MPI_COMM_WORLD);

            /* Orthogonalize first block using parallel CA-TSQR */
            R_ = malloc((s+1)*(s+1)*sizeof(double)); 
            TSQR_on_transpose(V, m, s, R_, myid, nprocs, MPI_COMM_WORLD); // note: resulting R is transposed!
            MPI_Barrier(MPI_COMM_WORLD);

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
                free(B_)
            }
            free(R_);
        }else{
            /* Matrix powers kernel (note: saved as transpose - vectors in rows!) */
            V = malloc(s*m*sizeof(double));
            matrix_powers(A, v, V, s, m, myid, nprocs, MPI_COMM_WORLD);

            /* Block-CGS to orthogonalize compared to previous blocks */
            mathcalR_ = malloc((block*s + 1)*s*sizeof(double));
            bgs_on_transpose(mathcalQ, V, mathcalR_, m, block*s + 1, s, MPI_COMM_WORLD);
            // Note: mathcal R is not transposed

            /* Orthogonalize block using parallel CA-TSQR */
            R_ = malloc(s*s*sizeof(double)); 
            TSQR_on_transpose(V, m, s, R_, myid, nprocs, MPI_COMM_WORLD); // note: resulting R is transposed!
            MPI_Barrier(MPI_COMM_WORLD); 

            /* Set mathcal Q (note: saved as transpose - vectors in rows!) */
            memcpy(mathcalQ + (1 + block*s)*m, V, s*m*sizeof(double));

            /* Save last vector in v */
            memcpy(v, V + (s-1)*m, m*sizeof(double));
            free(V);

            /* Update mathcal H */
            update_hess_on_transpose(&mathcalH, mathcalR_, R_, s, block);
        }
    }
    
    MPI_Finalize();
    return 0;
}