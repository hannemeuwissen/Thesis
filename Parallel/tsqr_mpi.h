/**
 * @file tsqr_mpi.h
 * @brief Header for thesis at Trinity College Dublin.
 * @author Hanne Meuwissen (22307813)
 * @version 2.0
 * @date 2023-05-27
 */
#ifndef TSQR_H_7BZWC1IU
#define TSQR_H_7BZWC1IU

void TSQR(double *A, const int m, const int N, double *R, const int rank, const int nprocs, MPI_Comm comm);

#endif /* end of include guard: TSQR_H_7BZWC1IU */