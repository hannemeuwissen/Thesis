/**
 * @file sparse.c
 * @brief Code for parallel calculations with sparse CSR matrices.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 1.0
 * @date 2023-06-02
 */
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<mpi.h>
#include"sparse.h"
#include<mkl_types.h>
#include<mkl_cblas.h>
#include"mkl.h"

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
 * @brief Get the start and end indices of the rows for each part that each process holds of
 * the full matrix. 
 * @param n The total number of rows.
 * @param nprocs The number of processes.
 * @param start The start indices.
 * @param end The end indices.
 */
void get_indices(const int n, const int nprocs, int * start, int * end){
    int n_elements = n/nprocs;
    int remainder = n%nprocs;
    start[0] = 0;
    end[0] = ((!remainder) ? (n_elements-1) : n_elements);
    if(nprocs > 1){
        for(int i=1;i<nprocs;i++){
            start[i] = end[i-1] + 1;
            end[i] = end[i-1] + ((i<remainder) ? (n_elements+1) : n_elements); 
        }
    }

}

/**
 * @brief Function that finds which rank has the vector value matching a given column index.
 * @param colindex The index of the to be found element.
 * @param nprocs The total number of processes.
 * @param M The total length of the vector.
 * @param smaller Indicator if the rank of the process where the elements is stored is smaller than the
 *                calling rank.
 * @param rank The rank of the calling process.
 * @return index_data The rank that stores the element and the index of the element.
 */
// index_data find_rank_colindex(const int colindex, const int nprocs, int * start, int * end, const int smaller, const int rank){
//     index_data result;
//     result.rank = -1;
//     result.index = -1;
//     int ll, ul;
//     if(smaller){
//         ll = 0;
//         ul = rank;
//     }else{
//         ll = rank+1;
//         ul = nprocs;
//     }
//     for(int i=ll;i<ul;i++){
//         if(colindex <= end[i]){
//             result.rank = i;
//             result.index = colindex - start[i];
//             return result;
//         }
//     }
//     return result;
// }
int find_rank_colindex(const int colindex, const int nprocs, int * end, const int smaller, const int rank){
    int result = -1;
    int ll, ul;
    if(smaller){
        ll = 0;
        ul = rank;
    }else{
        ll = rank+1;
        ul = nprocs;
    }
    for(int i=ll;i<ul;i++){
        if(colindex <= end[i]){
            return i;
        }
    }
    return result;
}

/**
 * @brief Function that calculates the sparse matrix mector multiplication between
 * a sparse_CSR matrix and a vector with compatible length.
 * @param A Sparse CSR matrix structure. Part of the matrix stored in memory of the calling process.
 * @param x Vector. Part of the full vector that is stored in memory of the calling process.
 * @param len Length of the vector part. 
 * @param result Resulting part of the full vector.
 * @param myid Id of the calling process.
 * @param nprocs Number of processes.
 * @param comm MPI communicator between processes.
 */
void spmv(sparse_CSR A, double * x, double len, double * result, const int myid, const int nprocs, MPI_Comm comm){
    
    if(len != A.nrows){
        perror("Incompatible dimensions in parallel spmv.\n");
        exit(EXIT_FAILURE);
    }
    A.
    int M = A.ncols;
    // double x_element;
    double * x_gathered_elements = malloc(M*sizeof(int)); // maximum size
    int * start = malloc(nprocs*sizeof(int));
    int * end = malloc(nprocs*sizeof(int));
    get_indices(M, nprocs, start, end);

    MPI_Win win;
    MPI_Win_create(x, len*sizeof(double), sizeof(double), MPI_INFO_NULL, comm, &win);
    
    for(int i=0;i<len;i++){
        // result[i] = 0.0;
        int nnz_i = 0;
        for(int j=A.rowptrs[i];j<A.rowptrs[i+1];j++){
            // x_element = 0.0;
            int colindex = A.colindex[j];

            MPI_Win_fence(MPI_MODE_NOPRECEDE | MPI_MODE_NOPUT | MPI_MODE_NOSTORE,win);

            if(colindex >= start[myid] && colindex <= end[myid]){ /* Element from x in own memory */
                // x_element = x[colindex - start[myid]];
                x_gathered_elements[nnz_i] = x[colindex - start[myid]];
            }else{ /* Element from x in other processes' memory*/
                int smaller = ((colindex < start[myid]) ? 1 : 0);
                int colindex_rank = find_rank_colindex(colindex, nprocs, end, smaller, myid);
                MPI_Get(x_gathered_elements + nnz_i, 1, MPI_DOUBLE, colindex_rank, colindex - start[colindex_rank], 1, MPI_DOUBLE, win);
            }

            MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE | MPI_MODE_NOPUT,win);

            // result[i] += A.values[j]*x_element;
            nnz_i++;
        }

        result[i] = cblas_ddot(nnz_i, A.values + A.rowptrs[i], 1, x_gathered_elements, 1);
    }

    MPI_Win_free(&win);

    free(start);
    free(end);
    free(x_gathered_elements);
}
// void spmv(sparse_CSR A, double * x, double len, double * result, const int myid, const int nprocs, MPI_Comm comm){
    
//     if(len != A.nrows){
//         perror("Incompatible dimensions in parallel spmv.\n");
//         exit(EXIT_FAILURE);
//     }
    
//     int M = A.ncols;
//     int * start = malloc(nprocs*sizeof(int));
//     int * end = malloc(nprocs*sizeof(int));
//     get_indices(M, nprocs, start, end);
//     double * x_gathered_elements = malloc(A.nnz*sizeof(double));

//     MPI_Datatype indexed_values;
//     MPI_Win win;
//     MPI_Win_create(x, len*sizeof(double), sizeof(double), MPI_INFO_NULL, comm, &win);
    
//     for(int i=0;i<len;i++){
//         int j = A.rowptrs[i];
//         int nnz_i = 0;
        
//         while(j<A.rowptrs[i+1]){
//             int colindex = A.colindex[j];
//             for(int p=0;p<nprocs;p++){ /* Per process: get all elements at once */
//                 int cnt = 0; /* number of elements to get from this process */
//                 int * indices = malloc((end[p] - start[p]+1)*sizeof(int)); /* Indices of element in other process */
//                 while(colindex < end[p] && j < A.rowptrs[i+1] && nnz_i < A.nnz){
//                     if(p != myid){
//                         indices[cnt] = colindex - start[p];
//                         cnt++;
//                     }else{
//                         x_gathered_elements[nnz_i] = x[colindex - start[p]];
//                     }
//                     nnz_i++;
//                     j++;
//                     if(j < A.rowptrs[i+1]){
//                         colindex = A.colindex[j];
//                     }
//                 }
//                 if(p != myid){
//                     int * lengths = malloc(cnt*sizeof(int));
//                     for(int i=0;i<cnt;i++){lengths[i] = 1;}
//                     MPI_Type_indexed(cnt, lengths, indices, MPI_DOUBLE, &indexed_values);
//                     MPI_Type_commit(&indexed_values);

//                     MPI_Win_fence(MPI_MODE_NOPRECEDE | MPI_MODE_NOPUT | MPI_MODE_NOSTORE,win);
//                     MPI_Get(x_gathered_elements + start[p], cnt, MPI_DOUBLE, p, 0, 1, indexed_values, win);
//                     MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE | MPI_MODE_NOPUT,win);

//                     free(lengths);
//                 }
//                 free(indices);
//             }
//             j++;

//         }

//         result[i] = cblas_ddot(A.nnz, A.values, 1, x_gathered_elements, 1);
//     }

//     MPI_Win_free(&win);
//     MPI_Type_free(&indexed_values);

//     free(start);
//     free(end);
//     free(x_gathered_elements);
// }
