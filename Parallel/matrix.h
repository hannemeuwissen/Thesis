/**
 * @file matrix.h
 * @brief Header file related to working with dense matrices, part of Thesis project in 
 * High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 2.0
 * @date 2023-06-02
 */

#ifndef MATRIX_H_TWYFUKH2
#define MATRIX_H_TWYFUKH2
double average(double * v, const int n);
void print_vector(double * v, const int len);
void print_matrix(double * A, const int n, const int m);
void save_matrix(double * A, const int n, const int m, char * filename);
void print_matrix_transposed(double * A, const int n, const int m);
void read_matrix_from_file(const char *const filename, const int skip, double *A, const int M, const int N);
void GatherQ(double * Q, const int m, const int n, const int myid, const int nprocs, int * start, MPI_Comm comm);
#endif /* end of include guard: MATRIX_H_TWYFUKH2 */
