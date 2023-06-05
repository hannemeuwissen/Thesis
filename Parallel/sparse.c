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
index_data find_rank_colindex(const int colindex, const int nprocs, const int M, const int smaller, const int rank){
    int start, end;
    index_data result;
    result.rank = -1;
    result.index = -1;
    if(smaller){
        for(int i=0;i<rank;i++){
            decomp1d(M, nprocs, i, &start, &end);
            if(colindex <= end){
                result.rank = i;
                result.index = colindex - start;
                return result;
            }
        }
    }else{
        for(int i=(rank + 1);i<nprocs;i++){
            decomp1d(M, nprocs, i, &start, &end);
            if(colindex <= end){
                result.rank = i;
                result.index = colindex - start;
                return result;
            }
        }
    }
    return result;
}

/**
 * @brief Function that calculates the sparse matrix mector multiplication between
 * a sparse_CSR matrix and a vector with compatible length.
 * @param A Sparse CSR matrix structure. Part of the matrix stored in memory of the calling process.
 * @param M Number of columns in the original problem.
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
    int start, end;
    decomp1d(M, nprocs, myid, &start, &end); /* Partition M rows over processes */

    MPI_Win win;
    MPI_Win_create(x, len*sizeof(double), sizeof(double), MPI_INFO_NULL, comm, &win);
    
    for(int i=0;i<len;i++){
        result[i] = 0.0;
        for(int j=A.rowptrs[i];j<A.rowptrs[i+1];j++){
            int colindex = A.colindex[j];
            double x_element;

            MPI_Win_fence(MPI_MODE_NOPUT,win);

            if(colindex >= start && colindex < start + len){ /* Element from x in own memory */
                x_element = x[colindex];
            }else{ /* Element from x in other processes' memory*/
                int smaller = ((colindex < start) ? 1 : 0);
                index_data colindex_data = find_rank_colindex(colindex, nprocs, M, smaller, myid);
                MPI_Get(&x_element, 1, MPI_DOUBLE, colindex_data.rank, colindex_data.index*sizeof(double), 1, MPI_DOUBLE, win);
            }
            MPI_Win_fence(MPI_MODE_NOSTORE,win);

            result[i] += A.values[j]*x_element;
            if(!myid){
                printf("Element: %lf\n", x_element);
            }
        }
    }

    MPI_Win_free(&win);
}
