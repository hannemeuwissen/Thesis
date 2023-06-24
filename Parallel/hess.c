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