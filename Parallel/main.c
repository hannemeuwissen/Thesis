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
#include"graph.h"
#include"sparse.h"
#include"tsqr_mpi.h"
#include"mkl.h"

int main(int argc, char **argv)
{  
    int myid, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(!myid){
        float logprocs = log2(nprocs);
        if(ceil(logprocs) != floor(logprocs)){
            fprintf(stderr,"Error: The number of processes needs to be a power of 2 (because of TSQR).\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    
    /* Test SPMV */
    int m = 10000;
    int n = 2500;
    int nnz_per_row = 2000;
    sparse_CSR M = generate_regular_graph_part_csr(n, m, nnz_per_row);
    printf("Process %d finished generating graph part of size %dx%d.\n", myid, n, m);
    // /* Print in order */
    // if(!myid){
    //     printf("Rank %d:\n", myid);
    //     print_CSR(&M);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(myid == 1){
    //     printf("Rank %d:\n", myid);
    //     print_CSR(&M);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(myid == 2){
    //     printf("Rank %d:\n", myid);
    //     print_CSR(&M);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(myid == 3){
    //     printf("Rank %d:\n", myid);
    //     print_CSR(&M);
    // }
    double * x = malloc(n*sizeof(double));
    for(int i=0;i<n;i++){x[i] = 1.0;}
    double * result = malloc(n*sizeof(double));
    double t1 = MPI_Wtime();
    spmv(M, x, n, result, myid, nprocs, MPI_COMM_WORLD);
    double t2 = MPI_Wtime();
    if(!myid){
        printf("First 50 lines from result on process 0:\n");
        print_vector(result, 50); // result should be 1 overall (sum of row elements)
        printf("Runtime: %lf\n", t2-t1);
    }

    // /* Test TSQR */
    // double A[20] = {9, 2, 7, 5, 10, 9, 8, 10, 9, 9, 8, 3, 3, 0, 4, 10, 6, 7, 10, 0};
    // double * R = malloc(4*4*sizeof(double));
    // TSQR(A, 5*nprocs, 4, R, myid, nprocs, MPI_COMM_WORLD);
    // if(!myid){
    //     print_matrix(R, 4, 4);
    // }
    // MPI_Bcast(R, 4*4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // cblas_dtrsm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 5, 4, 1.0, R, 4, A, 4);
    // /* Print in order */
    // if(!myid){
    //     printf("Rank %d:\n", myid);
    //     print_matrix(A, 5, 4);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(myid == 1){
    //     printf("Rank %d:\n", myid);
    //     print_matrix(A, 5, 4);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(myid == 2){
    //     printf("Rank %d:\n", myid);
    //     print_matrix(A, 5, 4);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(myid == 3){
    //     printf("Rank %d:\n", myid);
    //     print_matrix(A, 5, 4);
    // }

    // Read input: degree of Krylov subspace

    // Decide on s in s-step process

    // For all blocks:
    // 1. Matrix powers kernel using parallel spmv 
    // 2. Block-GS to orthogonalize compared to previous blocks (not the first time)
    // 3. Orthogonalize block using parallel CA-TSQR
    
    MPI_Finalize();
    return 0;
}