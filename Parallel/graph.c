/**
 * @file graph.c
 * @brief Code related to the generation of regular graph transition matrices in sparse 
 * format to investigate CA-Arnoldi, part of Thesis project in High Performance Computing 
 * at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 3.0
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
 * @brief Function that returns a result array with randomly selected and sorted column indices.
 * @param result Pointer to the address of the result array.
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
 * @brief Function that generates a part of the transition matrix (CSR) of a random graph 
 * where each node has the same number of edges, which are randomly selected.
 * @param n The number of nodes (rows) in the part.
 * @param M Total number of nodes in the graph.
 * @param nnz_per_row The number of nonzeros; the number of edges per node.
 * @param random Indicator if it generates using random device seeding or a fixed seed.
 * @return Sparse CSR matrix stucture of the part of the transition matrix.
 */
sparse_CSR generate_regular_graph_part_csr(const int n, const int M, const int nnz_per_row, const int random){
    if(!random){
        srandom(9499);
    }else{/* Seed random using random device */
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
    }

    /* Initialize sparse_CSR structure */
    sparse_CSR T;
    T.nrows = n;
    T.ncols = M;
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
        random_col_indices(&col_indices, M, nnz_per_row);
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

/**
 * @brief Function that generates a vector of random numbers of nonzeros between a minimim and maximum.
 * @param result The resulting vector.
 * @param sum The sum of all nonzeros.
 * @param min_nnz_per_row The maximum number of nonzeros per row.
 * @param max_nnz_per_row The minimum number of nonzeros per row.
 * @param n The number of rows.
 */
void random_nnz_per_row(int *result, int *sum, const int min_nnz_per_row, const int max_nnz_per_row, const int n){
    *sum = 0;
    for(int i=0;i<n;i++){
        result[i] = (min_nnz_per_row + random()) % (max_nnz_per_row+1 - min_nnz_per_row);
        *sum += result[i];
    }
}

/**
 * @brief Function that generates a part of the transition matrix (CSR) of a random graph 
 * where each node has a number of edges between specified min and max, which are randomly selected.
 * @param n The number of nodes (rows) in the part.
 * @param M Total number of nodes in the graph.
 * @param min_nnz_per_row The minimum of the number of nonzeros; the number of edges per node.
 * @param max_nnz_per_row The maximum of the number of nonzeros; the number of edges per node.
 * @param random_ind Indicator if it generates using random device seeding or a fixed seed.
 * @return Sparse CSR matrix stucture of the part of the transition matrix.
 */
sparse_CSR generate_irregular_graph_part_csr(const int n, const int M, const int min_nnz_per_row, const int max_nnz_per_row, const int random_ind){
    if(!random_ind){
        srandom(9499);
    }else{/* Seed random using random device */
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
    }

    /* Initialize sparse_CSR structure */
    sparse_CSR T;
    T.nrows = n;
    T.ncols = M;
    T.rowptrs = malloc((n+1)*sizeof(int));
    int *nnz_per_row = malloc(n*sizeof(int));
    random_nnz_per_row(nnz_per_row, &(T.nnz), min_nnz_per_row, max_nnz_per_row, n);
    T.colindex = malloc(T.nnz*sizeof(int));
    T.values = malloc(T.nnz*sizeof(double));

    /* Set nonzero elements per row at random */
    int i = 0;
    int row_index = 0;
    int * col_indices;
    double value;
    while(i<T.nnz){
        value = 1.0/((double) nnz_per_row[row_index]);
        T.rowptrs[row_index] = i;
        random_col_indices(&col_indices, M, nnz_per_row[row_index]);
        for(int j=0;j<nnz_per_row[row_index];j++){
            T.values[i] = value;            
            T.colindex[i] = col_indices[j];
            i++;
        }
        row_index++;
    }
    while(row_index <= n){ /* Account for zero nnz in last (few) rows */
       T.rowptrs[row_index] = T.nnz;
       row_index++;
    }
    free(col_indices);
    return T;
}

/**
 * @brief Function that generates a graph transition matrix (CSR) of a random graph 
 * where each node has a number of edges between a minimum and maximum, and the 
 * sparsity of the rows increases from the top to the bottom rows.
 * @param M Total number of nodes in the graph.
 * @param min_nnz The minimum of the number of nonzeros; the number of edges per node.
 * @param max_nnz The maximum of the number of nonzeros; the number of edges per node.
 * @return Sparse CSR matrix stucture of the part of the transition matrix.
 */
sparse_CSR generate_irregular_csr(const int M, const int min_nnz, const int max_nnz){
    srandom(9499);

    /* Initialize sparse_CSR structure */
    sparse_CSR T;
    T.nrows = M;
    T.ncols = M;
    T.rowptrs = malloc((M+1)*sizeof(int));
    int *nnz_per_row = malloc(M*sizeof(int));
    random_nnz_per_row(nnz_per_row, &(T.nnz), min_nnz, max_nnz, M);
    int sorter(const void * f1, const void * f2){return (*(int*)f2 - *(int*)f1);}
    qsort(nnz_per_row, M, sizeof(int), sorter);
    T.colindex = malloc(T.nnz*sizeof(int));
    T.values = malloc(T.nnz*sizeof(double));

    /* Set column indices nonzero elements per row at random */
    int i = 0;
    int row_index = 0;
    int * col_indices;
    double value;
    printf("last nr of nnz: %d\n", nnz_per_row[M-1]);
    while(i<T.nnz){
        value = 1.0/((double) nnz_per_row[row_index]);
        T.rowptrs[row_index] = i;
        random_col_indices(&col_indices, M, nnz_per_row[row_index]);
        for(int j=0;j<nnz_per_row[row_index];j++){
            T.values[i] = value;            
            T.colindex[i] = col_indices[j];
            i++;
        }
        row_index++;
    }
    while(row_index <= M){ /* Account for zero nnz in last (few) rows */
       T.rowptrs[row_index] = T.nnz;
       row_index++;
    }
    printf("Total nnz: %d\n", T.nnz);
    free(col_indices);
    free(nnz_per_row);
    return T;
}