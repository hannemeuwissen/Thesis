/**
 * @file matrix.h
 * @brief Header file for functions on dense matrices for thesis at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 1.0
 * @date 2023-06-02
 */

#ifndef MATRIX_H_TWYFUKH2
#define MATRIX_H_TWYFUKH2
void print_vector(double * v, const int len);
void print_matrix(double * A, const int n, const int m);
void read_matrix_from_file(const char *const filename, const int skip, double *A, const int M, const int N);
#endif /* end of include guard: MATRIX_H_TWYFUKH2 */
