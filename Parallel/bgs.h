/**
 * @file bgs.h
 * @brief Header file related to the calculation of block (classical) Gram-Schmidt 
 * as part of CA-Arnoldi, part of Thesis project in High Performance Computing at 
 * Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 2.0
 * @date 2023-06-11
 */
#ifndef BGS_H_7BZWC1IU
#define BGS_H_7BZWC1IU

void bgs(double *Q, double *V, double * R, const int m, const int n, const int s, MPI_Comm comm);
void bgs_on_transpose(double *Q, double *V, double * R, const int m, const int n, const int s, MPI_Comm comm);

#endif /* end of include guard: BGS_H_7BZWC1IU */