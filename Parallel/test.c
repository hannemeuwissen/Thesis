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
#include<string.h>
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

    // if(!myid){
    //     float logprocs = log2(nprocs);
    //     if(ceil(logprocs) != floor(logprocs)){
    //         fprintf(stderr,"Error: The number of processes needs to be a power of 2 (because of TSQR).\n");
    //         MPI_Abort(MPI_COMM_WORLD, 1);
    //     }
    // }

    // /* Test irregular graph generator */
    // sparse_CSR A = generate_irregular_graph_part_csr(5, 10, 2, 7, 1);
    // /* Print in order */
    // if(!myid){
    //     printf("Rank %d:\n", myid);
    //     print_CSR(&A);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(myid == 1){
    //     printf("Rank %d:\n", myid);
    //     print_CSR(&A);
    // }
    
    /* Test SPMV */
    int m = 20;
    int nnz_per_row = 5;
    int * start = malloc(nprocs*sizeof(int));
    int * end = malloc(nprocs*sizeof(int));
    get_indices(m, nprocs, start, end);
    int n = end[myid] - start[myid] + 1;
    sparse_CSR M = generate_regular_graph_part_csr(n, m, nnz_per_row, 1);
    // printf("Process %d finished generating graph part of size %dx%d.\n", myid, n, m);
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
    spmv(M, x, n, result, myid, nprocs, start, end, MPI_COMM_WORLD);
    double t2 = MPI_Wtime();
    printf("Process %d finished spmv\n",myid);
    if(!myid){
        printf("First lines from result on process 0:\n");
        print_vector(result, n); // result should be 1 overall (sum of row elements)
        printf("Runtime: %lf\n", t2-t1);
    }
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(myid == 3){
    //     printf("First lines from result on process 3:\n");
    //     print_vector(result, 10); // result should be 1 overall (sum of row elements)
    //     printf("Runtime: %lf\n", t2-t1);
    // }

    // /* Test TSQR */
    // int m = 25000; // Total: 100000
    // int n = 1000;
    // double * A = malloc(m*n*sizeof(double));
    // int skip = ((!myid) ? 0 : myid*m*n + n);
    // read_matrix_from_file("A.txt", skip, A, m, n); // Change skip: read_matrix_function changed 
    // if(myid == 3){
    //     printf("Data read!\n");
    // }
    // double * R = malloc(n*n*sizeof(double));
    // double t1 = MPI_Wtime();
    // TSQR(A, m, n, R, myid, nprocs, MPI_COMM_WORLD);
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
    // /* Test TSQR */
    // int m = 5; // Total: 4000
    // int n = 4;
    // double * A = malloc(m*n*sizeof(double));
    // double * transA = malloc(m*n*sizeof(double));
    // read_matrix_from_file("smallA.txt", 0, A, m, n); // Change skip: read_matrix_function changed 
    // // Transpose A
    // for(int i=0;i<m;i++){
    //     for(int j = 0;j<n;j++){
    //         transA[j*m + i] = A[i*n + j];
    //     }
    // }
    // // if(myid == 0){
    // //     print_matrix(transA, n, m);
    // // }
    // double * transR = malloc(n*n*sizeof(double));
    // double t1 = MPI_Wtime();
    // TSQR_on_transpose(transA, m, n, transR, myid, nprocs, MPI_COMM_WORLD);
    // double t2 = MPI_Wtime();
    // if(!myid){
    //     printf("Result for transR:\n");
    //     print_matrix(transR, n, n);
    //     printf("Result for transA:\n");
    //     print_matrix(transA, n, m);
    //     printf("Runtime: %lf\n", t2-t1);
    // }

    // /* Tesy BGS: 2 processes*/
    // double *V = malloc(6*sizeof(double)); 
    // double *W = malloc(4*sizeof(double));
    // memset(V, 0, 6*sizeof(double));
    // memset(W, 0, 4*sizeof(double));
    // if(!myid){
    //     V[0] = 1.0;
    //     V[4] = 1.0;
    //     W[0] = 1.0;
    //     W[3] = 5.0;
    // }else{
    //     V[2] = 1.0;
    //     V[4] = 1.0;
    //     W[0] = 3.0;
    //     W[1] = 7.0;
    //     W[2] = 4.0;
    // }
    // bgs(V, W, 2, 3, 2, MPI_COMM_WORLD);
    // if(!myid){
    //     print_matrix(W, 2, 2);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(myid == 1){
    //     print_matrix(W, 2, 2);
    // }   

    // /* Test read CSR */
    // sparse_CSR A;
    // read_CSR(&A, "smallcsr.txt");

    // print_CSR(&A); 

    MPI_Finalize();
    return 0;
}