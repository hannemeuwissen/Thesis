/**
 * @file graph.c
 * @brief Code for generating and doing calculations with graphs 
 * in sparse format.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 1.0
 * @date 2023-05-26
 */
#include<stdlib.h>
#include<stdio.h>
#include"graph.h"

int contains(size_t element, size_t * array, const size_t len){
    if(len > 0){
        for(int i=0;i<len;i++){
            if(array[i] == element){
                return 1;
            }
        }
    }
    return 0;
}

size_t * random_col_indices(const size_t nnz, const size_t n){
    srand(1999); // remove this once the code works!!!
    size_t * result = malloc(nnz*sizeof(size_t));
    size_t proposal;
    for(size_t i=0;i<nnz;i++){
        proposal = rand() % n;
        while(contains(proposal, result, i)){
            proposal = rand() % n;
        }
        result[i] = proposal;
    }
    return result;
}

sparse_CSR generate_regular_graph_trans_csr(const size_t n, const size_t nnz_per_row){

    /* Initialize sparse_CSR structure */
    sparse_CSR T;
    T.nrows = n;
    T.ncols = n;
    T.nnz = n*nnz_per_row;
    T.rowptrs = malloc((n+1)*sizeof(size_t));
    T.colindex = malloc(T.nnz*sizeof(size_t));
    T.values = malloc(T.nnz*sizeof(double));

    /* Set nonzero elements per row at random */
    size_t i = 0;
    size_t row_index = 0;
    double value = 1.0/((double) nnz_per_row)
    size_t * col_indices = random_col_indices(nnz_per_row);
    while(i<T.nnz){
        T.rowptrs[row_index++] = i;
        for(int j=0;j<nnz_per_row;j++){
            T.values[i] = value;            
            T.colindex[i] = col_indices[j];
            i++;
        }
        free(col_indices);
    }
}

void print_CSR(sparse_CSR * M){
    printf("Row pointers: ")
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