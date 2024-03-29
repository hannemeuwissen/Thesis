/**
 * @file graph.c
 * @brief Code related to generating sparse transition matrices of random graphs, 
 * part of Thesis project in High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 1.0
 * @date 2023-05-26
 */
#include<stdlib.h>
#include<stdio.h>
#include"graph.h"
#include"sparse.h"

/**
 * @brief Function that returns 1 if an index is allready present in a given array.
 * @param element Index.
 * @param array Array of previously proposed indices.
 * @param len Length of the array.
 * @return 1 if the array already contains the element, otherwise 0.
 */
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

/**
 * @brief Function that returns a result array with random selected and sorted column indices.
 * @param result The address of the result array.
 * @param n The number of columns.
 * @param nnz The number of nonzeros; the number of indices to be randomly selected.
 */
void random_col_indices(int ** result, const int n, const int nnz){
    *result = malloc(nnz*sizeof(int));
    int proposal;
    for(int i=0;i<nnz;i++){
        proposal = random() % n;
        while(contains(proposal, *result, i)){
            proposal = random() % n;
        }
        (*result)[i] = proposal;
    }
    int sorter(const void * f1, const void * f2){return (*(int*)f1 - *(int*)f2);}
    qsort(*result, nnz, sizeof(int), sorter);
}

/**
 * @brief Function that generates the transition matrix (CSR) of a random graph 
 * where each node has the same number of edges, which are randomly selected.
 * @param n The number of nodes.
 * @param nnz_per_row The number of nonzeros; the number of edges per node.
 * @return sparse_CSR Sparse CSR matrix stucture of the transition matrix.
 */
sparse_CSR generate_regular_graph_trans_csr(const int n, const int nnz_per_row){
    /* Seed random using random device */
    int randomvalue;
    FILE * fpointer;
    if((fpointer=fopen("/dev/random","r")) == NULL){
        perror("Error opening random device");
        exit(EXIT_FAILURE);
    }
    if(fread(&randomvalue,sizeof(int),1,fpointer) != 1){
        perror("Error reading from random device");
        exit(EXIT_FAILURE);
    }
    fclose(fpointer);
    srandom(randomvalue);

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
    int * col_indices;
    while(i<T.nnz){
        T.rowptrs[row_index++] = i;
        random_col_indices(&col_indices, n, nnz_per_row);
        for(int j=0;j<nnz_per_row;j++){
            T.values[i] = value;            
            T.colindex[i] = col_indices[j];
            i++;
        }
    }
    T.rowptrs[row_index] = T.nnz;
    free(col_indices);
    return T;
}