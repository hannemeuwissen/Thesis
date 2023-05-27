/**
 * @file sparse.c
 * @brief Code for calculations with sparse CSR matrices.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 1.0
 * @date 2023-05-27
 */
#include<stdlib.h>
#include<stdio.h>
#include"graph.h"

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