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

int contains(int element, int * array, const int len){
    if(len > 0){
        for(int i=0;i<len;i++){
            if(array[i] == element){
                return 1;
            }
        }
    }
    return 0;
}

int * random_col_indices(const int nnz, const int n){
    srand(1999); // remove this once the code works!!!
    int * result = malloc(nnz*sizeof(int));
    int proposal;
    for(int i=0;i<nnz;i++){
        proposal = rand() % n;
        while(contains(proposal, result, i)){
            proposal = rand() % n;
        }
        result[i] = proposal;
    }
    return result;
}

sparse_CSR generate_regular_graph_trans_csr(const int n, const int nnz_per_row){

    /* Initialize sparse_CSR structure */
    sparse_CSR T;
    T.nrows = n;
    T.ncols = n;
    T.nnz = n*nnz_per_row;
    T.rowptrs = malloc((n+1)*sizeof(int));
    T.colindex = malloc(T.nnz*sizeof(int));
    T.values = malloc(T.nnz*sizeof(double));

    /* Set nonzero elements per row at random */
    int i = 0;
    int row_index = 0;
    double value = 1.0/((double) nnz_per_row);
    int * col_indices = random_col_indices(nnz_per_row, n);
    while(i<T.nnz){
        T.rowptrs[row_index++] = i;
        for(int j=0;j<nnz_per_row;j++){
            T.values[i] = value;            
            T.colindex[i] = col_indices[j];
            i++;
        }
        // free(col_indices);
    }
    return T;
}

void print_CSR(sparse_CSR * M){
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