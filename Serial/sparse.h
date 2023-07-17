/**
 * @file sparse.h
 * @brief Header file related to working with sparse CSR matrices, part of Thesis 
 * project in High Performance Computing at Trinity College Dublin.
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

void read_CSR(sparse_CSR * M, const char *const filename);
void print_CSR(sparse_CSR * M);
void spmv(sparse_CSR M, double * v, double len, double * result);
void print_vector(double * v, const int len);
void print_matrix(double * A, const int n, const int m);
void read_matrix_from_file(const char *const filename, const int skip, double *A, const int M, const int N);
void Arnoldi(sparse_CSR A, double * b, const int len, double * Q, double * H, const int m);

#endif /* end of include guard: SPARSE_H_TWYFUKH2 */
