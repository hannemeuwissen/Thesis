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
#include"matrix.h"
#include"graph.h"
#include"sparse.h"
#include"tsqr_mpi.h"
#include"bgs.h"
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
    
    // /* Test SPMV */
    // int m = 10000;
    // int n = 2500;
    // int nnz_per_row = 2000;
    // sparse_CSR M = generate_regular_graph_part_csr(n, m, nnz_per_row);
    // printf("Process %d finished generating graph part of size %dx%d.\n", myid, n, m);
    // // /* Print in order */
    // // if(!myid){
    // //     printf("Rank %d:\n", myid);
    // //     print_CSR(&M);
    // // }
    // // MPI_Barrier(MPI_COMM_WORLD);
    // // if(myid == 1){
    // //     printf("Rank %d:\n", myid);
    // //     print_CSR(&M);
    // // }
    // // MPI_Barrier(MPI_COMM_WORLD);
    // // if(myid == 2){
    // //     printf("Rank %d:\n", myid);
    // //     print_CSR(&M);
    // // }
    // // MPI_Barrier(MPI_COMM_WORLD);
    // // if(myid == 3){
    // //     printf("Rank %d:\n", myid);
    // //     print_CSR(&M);
    // // }
    // double * x = malloc(n*sizeof(double));
    // for(int i=0;i<n;i++){x[i] = 1.0;}
    // double * result = malloc(n*sizeof(double));
    // double t1 = MPI_Wtime();
    // spmv(M, x, n, result, myid, nprocs, MPI_COMM_WORLD);
    // double t2 = MPI_Wtime();
    // if(!myid){
    //     printf("First 50 lines from result on process 0:\n");
    //     print_vector(result, 50); // result should be 1 overall (sum of row elements)
    //     printf("Runtime: %lf\n", t2-t1);
    // }

    // /* Test TSQR */
    // int m = 25000; // Total: 100000
    // int n = 1000;
    // double * A = malloc(m*n*sizeof(double));
    // int skip = ((!myid) ? 0 : myid*m*n + n);
    // read_matrix_from_file("A.txt", skip, A, m, n);
    // if(myid == 3){
    //     printf("Data read!\n");
    // }
    // double * R = malloc(n*n*sizeof(double));
    // double t1 = MPI_Wtime();
    // TSQR(A, m*nprocs, m, n, R, myid, nprocs, MPI_COMM_WORLD);
    // double t2 = MPI_Wtime();
    // if(!myid){
    //     printf("Result for R (first 10 elements of first row):\n");
    //     print_matrix(R, 1, 10);
    //     printf("Runtime: %lf\n", t2-t1);
    // }
    // t1 = MPI_Wtime();
    // MPI_Bcast(R, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // cblas_dtrsm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, 1.0, R, n, A, n);
    // MPI_Barrier(MPI_COMM_WORLD);
    // t2 = MPI_Wtime();
    // if(!myid){
    //     printf("Result for Q (first 10 elements of first row):\n");
    //     print_matrix(A, 1, 10);
    //     printf("Runtime: %lf\n", t2-t1);
    // }

    /* Tesy BGS: 2 processes*/
    double V[4], W[4];
    if(!myid){
        V = {1.0, 0.0, 0.0, 1.0};
        W = {1.0, 0.0, 0.0, 5.0};
    }else{
        V = {0.0, 0.0, 0.0, 1.0};
        W = {3.0, 7.0, 4.0, 0.0};
    }
    bgs(V, W, 4, 4, MPI_COMM_WORLD);
    if(!myid){
        print_matrix(W, 4, 4);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(myid == 1){
        print_matrix(W, 4, 4);
    }    

    MPI_Finalize();
    return 0;
}