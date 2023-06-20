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
    // char filename_A[100]; 
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
    // MPI_Bcast(filename_A, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(filename_v, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
    int steps = degree/s;
    int start, end;
    decomp1d(M, nprocs, myid, &start, &end);
    int m = end - start + 1;

    // Generate graph
    sparse_CSR A = generate_regular_graph_part_csr(m, M, nnz);

    // read v
    double *v = malloc(m*sizeof(double));
    read_matrix_from_file(filename_v, start, v, m, 1);

    // For all blocks:
    for(int block = 0;block < steps;block++){

        // Fix B

        // Matrix powers kernel using parallel spmv for size of block
        // NOTE: since output for each spmv is vector --> need to transpose before next part!
        int nr_vecs = ((!block)?s+1:s);
        double * V = malloc((s+1)*m*sizeof(double));
        if(!block){
            memcpy(V, v, m*sizeof(double)); // find way to always get the last vector from last block here
        }
        for(int k=0;k<s;k++){
            spmv(A, V + s*m, m, V + (s+1)*m, myid, nprocs, MPI_COMM_WORLD);
        }

        if(!block){
            // Orthogonalize block using parallel CA-TSQR
            double * R_ = malloc((s+1)*(s+1)*sizeof(double)); 
            TSQR(V, m, s, R_, myid, nprocs, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            // set mathcal Q
            double * mathcalQ = malloc(m*(s+1)*sizeof(double));
            memcpy(mathcalQ, V, m*(s+1)*sizeof(double)); // save mathcallQ in rows (transposed)?

            // Calculate mathcal H
            if(!myid){
                double * mathcalH = malloc((s+1)*s*sizeof(double));

                // set B_

                double * R = malloc(s*s*sizeof(double)); 
                get_R(R, R_, s+1);
                calc_hess(mathcalH, R_, B_, R, s+1, s);
            }
        }else{
            // Block-CGS to orthogonalize compared to previous blocks
            double mathcalR = malloc((block*s + 1)*s*sizeof(double));
            bgs(mathcalQ, V, mathcalR, m, block*s + 1, s, MPI_COMM_WORLD);

            // Orthogonalize block using parallel CA-TSQR
            double * R = malloc(s*s*sizeof(double)); 
            TSQR(V, m, s, R, myid, nprocs, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD); 

            // Update mathcal Q

            // Needed for following: mathcalR, mathcalR_, R
            // Compute Hk−1,k
            //Compute Hk viaEquation(3.21)
            // ComputehkviaEquation (3.22)
            // Hk = ... see p 165 Mark 

            // Update mathcal H
        }

    }
    
    MPI_Finalize();
    return 0;
}