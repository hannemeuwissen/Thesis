/**
 * @file bgs.c
 * @author Hanne Meuwissen (22307813)
 * @brief Code for thesis at Trinity College Dublin.
 * @version 0.1
 * @date 2023-06-11
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include<mkl_types.h>
#include<mkl_cblas.h>
#include<string.h>
#include"mkl.h"

/**
 * @brief Function that performs the block Gram-Schimdt update for W; orthogonalizes W w.r.t. V.
 * @param V Matrix V.
 * @param W Matrix W.
 * @param m The number of rows in the matrix part.
 * @param N The number of columns.
 * @param comm The MPI communicator.
 */
void bgs(double *V, double *W, const int m, const int N, MPI_Comm comm){
    /* Step 1: matrix-matrix multiplication */
    double * B = malloc(N*N*sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, N, m, 1.0, V, N, W, N, 0.0, B, N);

    /* Step 2: Communicate - Allreduce for each element */
    double * globalB = malloc(N*N*sizeof(double));
    MPI_Allreduce(B, globalB, N*N, MPI_DOUBLE, MPI_SUM, comm);
    free(B);

    /* Step 3: Local Gramm-Schmidt step - W <- W - VB */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, N, N, -1.0, V, N, globalB, N, 1.0, W, N);
    free(globalB);
}