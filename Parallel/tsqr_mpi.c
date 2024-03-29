/**
 * @file tsqr_mpi.c
 * @brief Code related to the calculation of CA-TSQR as part of CA-Arnoldi, part 
 * of Thesis project in High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 2.0
 * @date 2023-05-27
 */
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<mpi.h>
#include<mkl_types.h>
#include<mkl_cblas.h>
#include<string.h>
#include"mkl.h"
#include"sparse.h"
#include"matrix.h"

/**
 * @brief Function that determines if the process is active in a certain step of the TSQR algorithm.
 * @param rank Rank of the calling process.
 * @param step Step in the TSQR algorithm.
 * @return int Returns 0 if the calling process is not active anymore, otherwise 1.
 */
int is_active(const int rank, const int step){
    return (((int) rank%((int)pow(2.0, step))) == 0 ? 1 : 0);
}

/**
 * @brief Function that performs the communication avoiding TSQR (A = QR).
 * @param[in] A The matrix part of A for the calling process.
 * @param[in] m Number of rows in the matrix part.
 * @param[in] N Number of columns of the matrix.
 * @param[in] R The matrix that will hold the result for R. Rank 0 holds the final result.
 * @param[in] rank The rank of the calling process.
 * @param[in] nprocs The number of processes.
 * @param[in] comm The MPI communicator
 */
void TSQR(double *A, const int m, const int N, double *R, const int rank, const int nprocs, MPI_Comm comm){
    
    const int steps = log2(nprocs);
    
    double * tempA = malloc(m*N*sizeof(double)); /* Allocate space to not overwrite A */
    memcpy(tempA, A, m*N*sizeof(double));
    double * tau = malloc(N*sizeof(double)); /* Allocate space to hold tau's */

    for(int step=0;step<=steps;step++){ /* Parallel TSQR loop */

        if(is_active(rank, step)){

            lapack_int rows = ((!step)? m : 2*N); 
            lapack_int cols = N;

            /* Calculate local Housholder QR in terms of tau and v's */
            int ret = LAPACKE_dgeqrf(CblasRowMajor, rows, cols, tempA, cols, tau);
            if(ret!=0){
                if(ret<0){
                    fprintf(stderr, "LAPACKE_dgeqrf failed. Parameter %d had an illegal value\n", abs(ret));
                    exit(EXIT_FAILURE);
                }
                else{
                    fprintf(stderr, "LAPACKE_dgeqrf failed. Leading minor of order %d is not positive definite\n", ret);
                    exit(EXIT_FAILURE);
                }
            }

            /* Save R part */
            for(int i=0;i<N;i++){ 
                for(int j=0;j<N;j++){
                    R[j + i*N] = ((i>j) ? 0 : tempA[j + i*N]);
                }
            } 
            free(tempA);
            
            if(step<steps){
                if(is_active(rank, step + 1)){
                    /* Receive R from other process */
                    tempA = malloc(N*2*N*sizeof(double));
                    memcpy(tempA, R, N*N*sizeof(double));
                    MPI_Recv(tempA + N*N, N*N, MPI_DOUBLE, MPI_ANY_SOURCE, 1, comm, MPI_STATUS_IGNORE);
                }else{
                    /* Send R to other active process */
                    int lower_active = rank - pow(2, step);
                    MPI_Send(R, N*N, MPI_DOUBLE, lower_active, 1, comm);
                }
            }
        }else{
            break;
        }
    }
    free(tau);

    /* Broadcast result for R from process 0 */
    MPI_Bcast(R, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* Overwrite A with resulting Q for each part */
    cblas_dtrsm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, N, 1.0, R, N, A, N);
}

/**
 * @brief Function that performs the communication avoiding TSQR (A = QR), where the transpose of A is used.
 * @param[in] A The matrix part of A for the calling process (transposed).
 * @param[in] m Number of rows in the matrix part (before transposing).
 * @param[in] N Number of columns of the matrix (before transposing).
 * @param[in] R The matrix that will hold the result for R. Rank 0 holds the final result (transposed). 
 * @param[in] rank The rank of the calling process.
 * @param[in] nprocs The number of processes.
 * @param[in] comm The MPI communicator
 */
void TSQR_on_transpose(double *A, const int m, const int N, double *R, const int rank, const int nprocs, MPI_Comm comm){

    const int steps = log2(nprocs);
    
    double * tempA = malloc(m*N*sizeof(double)); /* Allocate space to not overwrite A */
    memcpy(tempA, A, m*N*sizeof(double));
    double * tau = malloc(N*sizeof(double)); /* Allocate space to hold tau's */

    for(int step=0;step<=steps;step++){ /* Parallel TSQR loop */

        if(is_active(rank, step)){

            lapack_int rows = ((!step)? m : 2*N); 
            lapack_int cols = N;

            /* Calculate local Housholder QR in terms of tau and v's */
            int ret = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, rows, cols, tempA, rows, tau);
            if(ret!=0){
                if(ret<0){
                    fprintf(stderr, "LAPACKE_dgeqrf failed. Parameter %d had an illegal value\n", abs(ret));
                    exit(EXIT_FAILURE);
                }
                else{
                    fprintf(stderr, "LAPACKE_dgeqrf failed. Leading minor of order %d is not positive definite\n", ret);
                    exit(EXIT_FAILURE);
                }
            }

            /* Save R part */
            for(int i=0;i<N;i++){ 
                for(int j=0;j<N;j++){
                    R[j + i*N] = ((j > i)? 0 : tempA[j + i*rows]);
                }
            }
            free(tempA);
            
            if(step<steps){
                if(is_active(rank, step + 1)){
                    /* Receive R from other process */
                    tempA = malloc(N*2*N*sizeof(double));
                    for(int i=0;i<N;i++){
                        memcpy(tempA + i*2*N, R + i*N, N*sizeof(double));
                    }
                    MPI_Recv(R, N*N, MPI_DOUBLE, MPI_ANY_SOURCE, 1, comm, MPI_STATUS_IGNORE);
                    for(int i=0;i<N;i++){
                        memcpy(tempA + N + i*2*N, R + i*N, N*sizeof(double));
                    }
                }else{
                    /* Send R to other active process */
                    int lower_active = rank - pow(2, step);
                    MPI_Send(R, N*N, MPI_DOUBLE, lower_active, 1, comm);
                }
            }
        }else{
            break;
        }
    }
    free(tau);

    /* Broadcast result for R from process 0 */
    MPI_Bcast(R, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* Overwrite A with resulting Q for each part */
    cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, N, 1.0, R, N, A, m);
}