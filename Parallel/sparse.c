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
#include"decomp1d.h"

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
 * @brief Function that finds which rank has the vector value matching a given column index.
 * @param colindex The index of the to be found element.
 * @param nprocs The total number of processes.
 * @param M The total length of the vector.
 * @param smaller Indicator if the rank of the process where the elements is stored is smaller than the
 *                calling rank.
 * @param rank The rank of the calling process.
 * @return index_data The rank that stores the element and the index of the element.
 */
// index_data find_rank_colindex(const int colindex, const int nprocs, const int M, const int smaller, const int rank){
//     int start, end;
//     index_data result;
//     result.rank = -1;
//     result.index = -1;
//     if(smaller){
//         for(int i=0;i<rank;i++){
//             decomp1d(M, nprocs, i, &start, &end);
//             if(colindex <= end){
//                 result.rank = i;
//                 result.index = colindex - start;
//                 return result;
//             }
//         }
//     }else{
//         for(int i=(rank + 1);i<nprocs;i++){
//             decomp1d(M, nprocs, i, &start, &end);
//             if(colindex <= end){
//                 result.rank = i;
//                 result.index = colindex - start;
//                 return result;
//             }
//         }
//     }
//     return result;
// }
index_data find_rank_colindex(const int colindex, const int nprocs, int * start, int * end, const int smaller, const int rank){
    index_data result;
    result.rank = -1;
    result.index = -1;
    if(smaller){
        for(int i=0;i<rank;i++){
            if(colindex <= end[i]){
                result.rank = i;
                result.index = colindex - start[i];
                return result;
            }
        }
    }else{
        for(int i=(rank + 1);i<nprocs;i++){
            if(colindex <= end[i]){
                result.rank = i;
                result.index = colindex - start[i];
                return result;
            }
        }
    }
    return result;
}

void get_indices(const int n, const int nprocs, int * start, int * end){
    int n_elements = n/nprocs;
    int remainder = n%nprocs;
    start[0] = 0;
    end[0] = ((!remainder) ? (n_elements) : (n_elements + 1));
    if(nprocs > 1){
        for(int i=1;i<nprocs;i++){
            start[i] = end[i-1] + 1;
            end[i] = end[i-1] + ((i<remainder) ? n_elements : (n_elements-1)); 
        }
    }

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
    
    int M = A.ncols;
    // int start, end;
    // decomp1d(M, nprocs, myid, &start, &end); /* Partition M rows over processes */
    double x_element;
    int * start = malloc(nprocs*sizeof(int));
    int * end = malloc(nprocs*sizeof(int));
    get_indices(M, nprocs, start, end);
    if(!myid){
        printf("Start: %d, end: %d\n", start[myid], end[myid]);
    }

    MPI_Win win;
    MPI_Win_create(x, len*sizeof(double), sizeof(double), MPI_INFO_NULL, comm, &win);
    
    for(int i=0;i<len;i++){
        result[i] = 0.0;
        for(int j=A.rowptrs[i];j<A.rowptrs[i+1];j++){
            x_element = 0.0;
            int colindex = A.colindex[j];

            MPI_Win_fence(MPI_MODE_NOPRECEDE | MPI_MODE_NOPUT | MPI_MODE_NOSTORE,win);

            if(colindex >= start[myid] && colindex <= end[myid]){ /* Element from x in own memory */
                x_element = x[colindex];
            }else{ /* Element from x in other processes' memory*/
                int smaller = ((colindex < start[myid]) ? 1 : 0);
                index_data colindex_data = find_rank_colindex(colindex, nprocs, start, end, smaller, myid);
                MPI_Get(&x_element, 1, MPI_DOUBLE, colindex_data.rank, colindex_data.index, 1, MPI_DOUBLE, win);
            }
            MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE | MPI_MODE_NOPUT,win);

            result[i] += A.values[j]*x_element;
        }
    }

    MPI_Win_free(&win);
}
