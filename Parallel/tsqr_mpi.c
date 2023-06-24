/**
 * @file tsqr_mpi.c
 * @author Hanne Meuwissen (22307813)
 * @brief Code for thesis at Trinity College Dublin.
 * @version 0.1
 * @date 2023-05-27
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include<mkl_types.h>
#include<mkl_cblas.h>
#include<string.h>
#include"mkl.h"
#include "sparse.h"
#include "matrix.h"

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
 * @brief Function that finds active process to send data to. 
 * @param rank Rank of calling process.
 * @param step Step in the TSQR algorithm.
 * @return int Returns the rank of the process it has to send its data to.
 */
int find_lower_active(const int rank, const int step){
    for(int i=rank-1;i>=0;i--){
        if(is_active(i, step)){
            return i;
        }
    }
    return 0;
}

/**
 * @brief Function that performs the communication avoiding TSQR (A = QR) using OpenMP.
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
                    int lower_active = find_lower_active(rank, step + 1);
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
    // MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * @brief Function that performs the communication avoiding TSQR (A = QR) using OpenMP.
 * @param[in] A The matrix part of A for the calling process.
 * @param[in] m Number of rows in the matrix part.
 * @param[in] N Number of columns of the matrix.
 * @param[in] R The matrix that will hold the result for R. Rank 0 holds the final result.
 * @param[in] rank The rank of the calling process.
 * @param[in] nprocs The number of processes.
 * @param[in] comm The MPI communicator
 */
void TSQR_on_transpose(double *A, const int m, const int N, double *R, const int rank, const int nprocs, MPI_Comm comm){
    
    MPI_Datatype stridedcol;
    MPI_Type_vector(N, N, m, MPI_DOUBLE, &stridedcol);
    MPI_Type_commit(&stridedcol);

    const int steps = log2(nprocs);
    
    double * tempA = malloc(m*N*sizeof(double)); /* Allocate space to not overwrite A */
    memcpy(tempA, A, m*N*sizeof(double));
    double * tau = malloc(N*sizeof(double)); /* Allocate space to hold tau's */

    for(int step=0;step<=steps;step++){ /* Parallel TSQR loop */

        if(is_active(rank, step)){

            lapack_int rows = ((!step)? m : 2*N); 
            lapack_int cols = N;

            /* Calculate local Housholder QR in terms of tau and v's */
            int ret = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, cols, rows, tempA, rows, tau);
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

            if(!rank){
                print_matrix(tempA, cols, rows);
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
                    MPI_Recv(tempA + N, N*N, MPI_DOUBLE, MPI_ANY_SOURCE, 1, comm, MPI_STATUS_IGNORE);
                }else{
                    /* Send R to other active process */
                    int lower_active = find_lower_active(rank, step + 1);
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
    cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, N, 1.0, R, N, A, N);
    // MPI_Barrier(MPI_COMM_WORLD);

    MPI_Type_free(&stridedcol);
}