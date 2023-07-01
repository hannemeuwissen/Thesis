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
 * @brief Function that performs the block Gram-Schimdt update for V; orthogonalizes V w.r.t. Q.
 * @param Q Matrix Q.
 * @param V Matrix V.
 * @param R Matrix R.
 * @param m The number of rows in the matrix part.
 * @param n The number of columns in Q.
 * @param s The number of columns in V
 * @param comm The MPI communicator.
 */
void bgs(double *Q, double *V, double * R, const int m, const int n, const int s, MPI_Comm comm){
    /* Step 1: matrix-matrix multiplication */
    double * localR = malloc(n*s*sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, s, m, 1.0, Q, n, V, s, 0.0, localR, s);

    /* Step 2: Communicate - Allreduce for each element */
    // double * globalR = malloc(n*s*sizeof(double));
    MPI_Allreduce(localR, R, n*s, MPI_DOUBLE, MPI_SUM, comm);
    free(localR);

    /* Step 3: Local Gramm-Schmidt step - V <- V - QR */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, s, n, -1.0, Q, n, R, s, 1.0, V, s);
    // free(globalR);
}

/**
 * @brief Function that performs the block Gram-Schimdt update for V; orthogonalizes V w.r.t. Q.
 * @param Q Matrix Q.
 * @param V Matrix V.
 * @param R Matrix R.
 * @param m The number of rows in the matrix part.
 * @param n The number of vectors in Q.
 * @param s The number of vectors in V
 * @param comm The MPI communicator.
 */
void bgs_on_transpose(double *Q, double *V, double * R, const int m, const int n, const int s, MPI_Comm comm){
    /* Step 1: matrix-matrix multiplication */
    double * localR = malloc(n*s*sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, s, m, 1.0, Q, n, V, s, 0.0, localR, n);

    /* Step 2: Communicate - Allreduce for each element */
    MPI_Allreduce(localR, R, n*s, MPI_DOUBLE, MPI_SUM, comm);
    free(localR);

    /* Step 3: Local Gramm-Schmidt step - V <- V - QR */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, s, n, -1.0, Q, m, R, s, 1.0, V, m);
}