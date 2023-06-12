/**
 * @file bgs.c
 * @author Hanne Meuwissen (22307813)
 * @brief Code for thesis at Trinity College Dublin.
 * @version 0.1
 * @date 2023-06-11
 */
#include <stdlib.h>
#include <stdio.h>
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
 * @param n The number of columns in V.
 * @param s The number of columns in W
 * @param comm The MPI communicator.
 */
void bgs(double *V, double *W, const int m, const int n, const int s, MPI_Comm comm){
    /* Step 1: matrix-matrix multiplication */
    double * B = malloc(n*s*sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, s, m, 1.0, V, n, W, s, 0.0, B, s);

    /* Step 2: Communicate - Allreduce for each element */
    double * globalB = malloc(n*s*sizeof(double));
    MPI_Allreduce(B, globalB, n*s, MPI_DOUBLE, MPI_SUM, comm);
    free(B);

    /* Step 3: Local Gramm-Schmidt step - W <- W - VB */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, s, n, -1.0, V, n, globalB, s, 1.0, W, s);
    free(globalB);
}