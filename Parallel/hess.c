/**
 * @file hess.c
 * @author Hanne Meuwissen (22307813)
 * @brief Code for thesis at Trinity College Dublin.
 * @version 0.1
 * @date 2023-06-19
 */
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include<mkl_types.h>
#include<mkl_cblas.h>
#include<string.h>
#include"mkl.h"

/**
 * @brief Function that calculates the upper Hessenberg matrix based on R_, R and B_. 
 * @param H_ The matrix which will hold the result.
 * @param R_ The R_ matrix.
 * @param B_ The B_ matrix.
 * @param R The R matrix.
 * @param n The number of rows and columns in R_, and rows in B_.
 * @param m The number of columns in B_, and rows and columns in R.
 */
void calc_hess(double * H_, double * R_, double * B_, double * R, const int n, const int m){
    memcpy(H_, B_, n*m*sizeof(double));
    cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, m, 1.0, R_, n, H_, m);
    cblas_dtrsm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n, m, 1.0, R, m, H_, m);
}

/**
 * @brief Function that calculates the upper Hessenberg matrix based on R_, R (both transposed) and B_. 
 * @param H_ The matrix which will hold the result.
 * @param R_ The R_ matrix.
 * @param B_ The B_ matrix.
 * @param R The R matrix.
 * @param n The number of rows and columns in R_, and rows in B_.
 * @param m The number of columns in B_, and rows and columns in R.
 */
void calc_hess_on_transpose(double * H_, double * R_, double * B_, double * R, const int n, const int m){
    memcpy(H_, B_, n*m*sizeof(double));
    cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, n, m, 1.0, R_, n, H_, m);
    cblas_dtrsm(CblasRowMajor, CblasRight, CblasUpper, CblasTrans, CblasNonUnit, n, m, 1.0, R, m, H_, m);
}

void update_hess_on_transpose(double ** H, double * mathcalR_, double * R, const int s, const int k){
    /* Set dimensions */
    const int lowerdim = s*(k+1);
    const int upperdim = s*(k+1) + 1;

    /* Construct MathcalR_ and MathcalR */
    double * MathcalR_ = malloc(upperdim*upperdim*sizeof(double));
    for(int i=0;i<(s*k + 1);i++){
        memset(MathcalR_ + i*upperdim, 0, (s*k + 1)*sizeof(double));
        MathcalR_[i + i*upperdim] = 1.0;
        memcpy(MathcalR_ + (s*k + 1) + i*upperdim, mathcalR_ + i*s, s*sizeof(double));
    }
    double * transR_ = malloc(s*s*sizeof(double));
    transpose(R_, transR_, s);
    for(int i=(s*k + 1);i<upperdim;i++){
        memset(MathcalR_ + i*upperdim, 0, (s*k + 1)*sizeof(double));
        memcpy(MathcalR_ + (s*k + 1) + i*upperdim, transR_ + i*s, s*sizeof(double));
    }
    double * MathcalR = malloc((s*k + 1)*(s*k + 1)*sizeof(double));
    get_R(MathcalR, MathcalR_, upperdim);

    /* Construct MathcalB_ */
    double * MathcalB_ = malloc(upperdim*lowerdim*sizeof(double));
    for(int i=0;i<(s*k + 1);i++){
        memcpy(MathcalB_ + i*lowerdim, *H + i*s*k, s*k*sizeof(double));
        memset(MathcalB_ + s*k + i*lowerdim, 0, s*sizeof(double));
    }
    /* Change of basis matrix */
    double * B_ = malloc(s*(s+1)*sizeof(double));
    double * B = malloc(s*s*sizeof(double));
    set_B(B, B_, s);

    /* Update H */
    double *tempH = malloc(upperdim*lowerdim*sizeof(double));
}

/**
 * @brief Function that gets the internal part from R_. 
 * @param R The matrix which will hold the result.
 * @param R_ The matrix which holds R_.
 * @param n The number of rows in R_.
 */
void get_R(double * R, double * R_, const int n){
    for(int i=0;i<(n-1);i++){
        memcpy(R + (n-1)*i, R_ + n*i, (n-1)*sizeof(double));
    }
}

/**
 * @brief Function that sets both change of basis matrices B_ and the internal part B. 
 * @param B The matrix which will hold the internal part of B_.
 * @param B_ The matrix which will hold B_.
 * @param s The number of columns in B_.
 */
void set_B(double *B, double *B_, const int s){
    /* Set B_ */    
    memset(B_, 0, s*(s+1)*sizeof(double));
    for(int i=1;i<(s+1);i++){
        for(int j=0;j<s;j++){
            if(i == (j+1)){
                B_[i*s + j] = 1.0;
            }
        }
    }

    /* Set B */
    memcpy(B, B_, s*s*sizeof(double));
}

/**
 * @brief Function that transposes a square matrix.
 * @param R Pointer to the matrix.
 * @param transR Poiter to matrix that stores the result.
 * @param s The dimension of the matrix.
 */
void transpose(double *R, double *transR, const int s){
    for(int i=0;i<s;i++){
        for(int j=0;j<s;j++){
            transR[j + i*s] = R[i + j*s];
        }
    }
}