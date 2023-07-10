/**
 * @file sparse.h
 * @brief Header file related to working with sparse CSR matrices, part of Thesis 
 * project in High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 2.0
 * @date 2023-05-26
 */
#ifndef SPARSE_H_TWYFUKH2
#define SPARSE_H_TWYFUKH2
#include<mpi.h>

/**
 * @brief Structure that holds information of sparse CSR matrix.
 */
typedef struct sparse_CSR {
    int nrows;
    int ncols;
    int nnz;
    int * rowptrs;
    int * colindex;
    double * values;
} sparse_CSR;

void print_CSR(sparse_CSR * M);
void read_CSR(sparse_CSR * M, const char *const filename);
void get_indices(const int n, const int nprocs, int * start, int * end);
void spmv(sparse_CSR A, double * x, double len, double * result, const int myid, const int nprocs, int * start, int * end, MPI_Comm comm, MPI_Win win);
void matrix_powers(sparse_CSR A, double * start_v, double * V, const int s, const int m, const int myid, const int nprocs, int *start, int *end, MPI_Comm comm);

#endif /* end of include guard: SPARSE_H_TWYFUKH2 */
