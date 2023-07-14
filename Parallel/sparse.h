/**
 * @file sparse.h
 * @brief Header file related to working with sparse CSR matrices, part of Thesis 
 * project in High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 5.0
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
void read_CSR_data(sparse_CSR * M, const char * filename);
void save_CSR(char * filename_A, sparse_CSR A);
void read_CSR_values(sparse_CSR * M, const char * filename, const int start, const int end)
void get_indices(const int n, const int nprocs, int * start, int * end);
void get_indices_load_balanced(sparse_CSR A, const int nprocs, int * start, int * end);
void spmv(sparse_CSR A, double * x, double len, double * result, const int myid, const int nprocs, int * start, int * end, MPI_Comm comm);
void matrix_powers(sparse_CSR A, double * start_v, double * V, const int s, const int m, const int myid, const int nprocs, int *start, int *end, MPI_Comm comm);

#endif /* end of include guard: SPARSE_H_TWYFUKH2 */
