/**
 * @file sparse.h
 * @brief Main function for thesis.
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

void print_CSR(sparse_CSR * M);
void spmv(sparse_CSR M, double * v, double len, double * result);
void print_vector(double * v, const int len);
void print_matrix(double * A, const int n, const int m)
void Arnoldi(sparse_CSR A, double * b, const int len, double * Q, const int m);

#endif /* end of include guard: SPARSE_H_TWYFUKH2 */
