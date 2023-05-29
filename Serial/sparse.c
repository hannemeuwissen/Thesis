/**
 * @file sparse.c
 * @brief Code for calculations with sparse CSR matrices.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 1.0
 * @date 2023-05-27
 */
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<mkl_types.h>
#include<mkl_cblas.h>
#include"sparse.h"

/**
 * @brief Function that prints a sparse CSR matrix structure.
 * @param M The sparse_CSR structure to be printed.
 */
void print_CSR(sparse_CSR * M){
    printf("Number of columns and rows: %d\nNumber of nonzeros: %d\n", M->ncols, M->nnz);
    printf("Row pointers: ");
    for(int i=0;i<(M->nrows+1);i++){
        printf("%d ", M->rowptrs[i]);
    }
    printf("\nColumn indices: ");
    for(int i=0;i<(M->nnz);i++){
        printf("%d ", M->colindex[i]);
    }
    printf("\nValues: ");
    for(int i=0;i<(M->nnz);i++){
        printf("%lf ", M->values[i]);
    }
    printf("\n");
}

/**
 * @brief Function that calculates the sparse matrix mector multiplication between
 * a sparse_CSR matrix and a vector with compatible length.
 * @param M Sparse CSR matrix structure.
 * @param v Vector.
 * @param len Length of the vector.
 * @param result Result vector.
 */
void spmv(sparse_CSR M, double * v, double len, double * result){
    if(len != M.nrows){
        perror("incompatible dimensions in spmv.\n");
        exit(EXIT_FAILURE);
    }
    for(int i=0;i<len;i++){
        result[i] = 0.0;
        for(int j=M.rowptrs[i];j<M.rowptrs[i+1];j++){
            result[i] += M.values[j]*v[M.colindex[j]];
        }
    }
}

/**
 * @brief Function that prints a vector of given length.
 * @param v Vector.
 * @param len Length of the vector.
 */
void print_vector(double * v, const int len){
    for(int i = 0;i<len;i++){
        printf("%lf\n", v[i]);
    }
}

/**
 * @brief Function that prints a matrix of given dimensions.
 * @param A The matrix.
 * @param n The number of rows.
 * @param m The number of columns.
 */
void print_matrix(double * A, const int n, const int m){
    for(int i = 0;i<n;i++){
        for(int j=0;j<m;j++){
            printf("%lf ", A[i*m + j]);
        }
        printf("\n");
    }
}

/**
 * @brief Function that computes the m-Krylov subspace of the sparse CSR matrix A,
 * thus the space spanned by the vectors {b, Ab, ..., A^(m-1) b}.
 * @param A Sparse sparse_CSR matrix struct. A should represent a len x len matrix.
 * @param b Vector b.
 * @param len Length of the vector.
 * @param Q Matrix which will hold the resulting basis vectors of the Krylov subspace. (len x m)
 * @param m The degree of the Krlov subspace.
 */
void Arnoldi(sparse_CSR A, double * b, const int len, double * Q, const int m){
    if(A.ncols != len){
        perror("Incompatible dimensions in Arnoldi.");
        exit(EXIT_FAILURE);
    }
    double eps = 1e-12;
    double * Q_trans = malloc(len*m*sizeof(double));
    cblas_dcopy(len, b, 1, Q_trans, 1); /* Set Q_trans0 */
    double norm_b = cblas_dnrm2(len, Q_trans, 1);
    cblas_dscal(len, 1/norm_b, Q_trans, 1); /* Normalize */
    double h;
    double * w = malloc(len*sizeof(double));
    for(int j=1;j < m;j++){
        spmv(A, Q_trans + (j-1)*len, len, w);
        print_vector(w, len);
        for(int i=0;i<j;i++){
            h = cblas_ddot(len, w, 1, Q_trans + i*len, 1);
            cblas_daxpy(len,-h,Q_trans + i*len,1,w, 1);
        }
        h = cblas_dnrm2(len, w, 1);
        if(h < eps){
            break;
        }
        cblas_dscal(len, 1/h, w, 1);
        cblas_dcopy(len, w, 1, Q_trans + j*len, 1);
    }
    /* Transpose final result */
    for(int i=0;i<len;i++){
        for(int j=0;j<m;j++){
            Q[i*m + j] = Q_trans[j*len + i];
        }
    }
}