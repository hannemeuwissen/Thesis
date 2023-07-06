/**
 * @file hess.c
 * @author Hanne Meuwissen (22307813)
 * @brief Code for thesis at Trinity College Dublin.
 * @version 0.1
 * @date 2023-06-19
 */
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<mpi.h>
#include<mkl_types.h>
#include<mkl_cblas.h>
#include<string.h>
#include"mkl.h"
#include"matrix.h"

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
    cblas_dtrmm(CblasRowMajor, CblasLeft, CblasLower, CblasTrans, CblasNonUnit, n, m, 1.0, R_, n, H_, m);
    cblas_dtrsm(CblasRowMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, n, m, 1.0, R, m, H_, m);
}

/**
 * @brief Function that sets the change of basis matrix B_. 
 * @param B_ The pointer to the matrix B_.
 * @param s The number of columns in B_.
 */
void set_B_(double *B_, const int s){
    memset(B_, 0, s*(s+1)*sizeof(double));
    for(int i=1;i<(s+1);i++){
        B_[i-1 + i*s] = 1.0;
    }
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

/**
 * @brief Function that updates the upper Hessenberg matrix from step 1.
 * @param H Pointer to the address of the previous upper Hessenberg matrix, and will store the next.
 * @param mathcalR_ Pointer to the matrix inproduct from B-CGS.
 * @param R Pointer to the R matrix from TSQR.
 * @param s The blocksize.
 * @param k The index of the current block.
 */
void update_hess_on_transpose(double ** H, double * mathcalR_, double * R_, const int s, const int k){
    /* Set dimensions */
    const int lowerdim = s*(k+1);
    const int upperdim = s*(k+1) + 1;

    /* Construct MathcalR_ ([I mathcalR_ ; 0 R_]) and MathcalR */
    double * MathcalR_ = malloc(upperdim*upperdim*sizeof(double));
    for(int i=0;i<(s*k + 1);i++){ /* Upper part: [I mathcalR_]*/
        memset(MathcalR_ + i*upperdim, 0, (s*k + 1)*sizeof(double));
        MathcalR_[i + i*upperdim] = 1.0;
        memcpy(MathcalR_ + (s*k + 1) + i*upperdim, mathcalR_ + i*s, s*sizeof(double));
    }
    double * transR_ = malloc(s*s*sizeof(double));
    transpose(R_, transR_, s);
    int j = 0;
    for(int i=(s*k + 1);i<upperdim;i++){ /* Lower part: [0 R_]*/
        memset(MathcalR_ + i*upperdim, 0, (s*k + 1)*sizeof(double));
        memcpy(MathcalR_ + (s*k + 1) + i*upperdim, transR_ + j*s, s*sizeof(double));
        j++;
    }
    double * MathcalR = malloc(lowerdim*lowerdim*sizeof(double));
    get_R(MathcalR, MathcalR_, upperdim);

    /* Construct MathcalB_ ([H 0 ; 0,...,h B_]) */
    double * MathcalB_ = malloc(upperdim*lowerdim*sizeof(double));
    for(int i=0;i<s*k;i++){ /* Upper part: [H 0] */
        memcpy(MathcalB_ + i*lowerdim, *H + i*s*k, s*k*sizeof(double));
        memset(MathcalB_ + s*k + i*lowerdim, 0, s*sizeof(double));
    }
    for(int i=s*k;i<upperdim;i++){ /* Lower part: [0,...,h B_] */
        memset(MathcalB_ + i*lowerdim, 0, lowerdim*sizeof(double));
        if(i > s*k){MathcalB_[i-1 + i*lowerdim] = 1.0;}
    }
    MathcalB_[s*k - 1 + s*k*lowerdim] = (*H)[s*k - 1 + (s*k)*(s*k)];
    free(*H);

    /* Update H */
    *H = malloc(upperdim*lowerdim*sizeof(double));
    calc_hess(*H, MathcalR_, MathcalB_, MathcalR, upperdim, lowerdim);

    free(MathcalB_);
    free(MathcalR_);
    free(MathcalR);
}

/**
 * @brief Function that checks if breakdown has occured in calculated Hessenberg
 * @param H The Hessenberg matrix.
 * @param s The blocksize.
 * @param k The current block index.
 * @return The function returns -1 when no breakdown has occured, and the row index of 
 * the first zero if breakdown occured.
 */
int breakdown_check(double *H, const int s, const int k, const double tol){
    printf("Breakdown check:\n")
    const int lowerdim = s*(k+1);
    const int upperdim = s*(k+1) + 1;

    for(int j=s*k;j<lowerdim;j++){
        printf("%lf\n", fabs(H[j + (j+1)*lowerdim]));
        if(fabs(H[j + (j+1)*lowerdim]) < tol){
            return j+1;
        }
    }
    return -1;
}