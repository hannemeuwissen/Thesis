/**
 * @file sparse.h
 * @brief Header file for thesis.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 1.0
 * @date 2023-05-26
 */
#ifndef SPARSE_H_TWYFUKH2
#define SPARSE_H_TWYFUKH2

typedef struct sparse_CSR {
    int nrows;
    int ncols;
    int nnz;
    int * rowptrs;
    int * colindex;
    double * values;
} sparse_CSR;

typedef struct index_data {
    int rank;
    int index;
} index_data;

void print_CSR(sparse_CSR * M);
void print_vector(double * v, const int len);
void print_matrix(double * A, const int n, const int m);
void spmv(sparse_CSR M, double * v, double len, double * result);

#endif /* end of include guard: SPARSE_H_TWYFUKH2 */
