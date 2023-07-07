/**
 * @file bgs.c
 * @brief Code related to the calculation of block (classical) Gram-Schmidt as part of CA-Arnoldi,
 * part of Thesis project in High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 2.0
 * @date 2023-06-11
 */
#include<stdlib.h>
#include<stdio.h>
#include<mpi.h>
#include<mkl_types.h>
#include<mkl_cblas.h>
#include<string.h>
#include"mkl.h"

/**
 * @brief Function that performs the block Gram-Schimdt update for V; orthogonalizes V w.r.t. Q.
 * @param Q Matrix Q.
 * @param V Matrix V.
 * @param R Matrix R, will hold the product Q'V.
 * @param m The number of rows in the matrix part (rows in Q).
 * @param n The number of columns in Q.
 * @param s The number of columns in V.
 * @param comm The MPI communicator.
 */
void bgs(double *Q, double *V, double * R, const int m, const int n, const int s, MPI_Comm comm){
    /* Step 1: matrix-matrix multiplication */
    double * localR = malloc(n*s*sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, s, m, 1.0, Q, n, V, s, 0.0, localR, s);

    /* Step 2: Communicate - Allreduce for each element */
    MPI_Allreduce(localR, R, n*s, MPI_DOUBLE, MPI_SUM, comm);
    free(localR);

    /* Step 3: Local Gramm-Schmidt step - V <- V - QR */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, s, n, -1.0, Q, n, R, s, 1.0, V, s);
}

/**
 * @brief Function that performs the block Gram-Schimdt update for V using the transpose of V and Q;
 * orthogonalizes V w.r.t. Q.
 * @param Q Matrix Q.
 * @param V Matrix V.
 * @param R Matrix R, will hold the product QV'.
 * @param m The number of rows in the matrix part (number of columns in Q).
 * @param n The number of vectors in Q (number of rows in Q).
 * @param s The number of vectors in V (number of rows in V).
 * @param comm The MPI communicator.
 */
void bgs_on_transpose(double *Q, double *V, double * R, const int m, const int n, const int s, MPI_Comm comm){
    /* Step 1: matrix-matrix multiplication */
    double * localR = malloc(n*s*sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, s, m, 1.0, Q, m, V, m, 0.0, localR, s);

    /* Step 2: Communicate - Allreduce for each element */
    MPI_Allreduce(localR, R, n*s, MPI_DOUBLE, MPI_SUM, comm);
    free(localR);

    /* Step 3: Local Gramm-Schmidt step - V <- V - QR */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, s, n, -1.0, Q, m, R, s, 1.0, V, m);
}