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
#include "decomp1d.h"
#include "sparse.h"

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
 * @brief Function that performs the communication avoiding TSQR using OpenMP.
 * @param[in] A The matrix that will hold the data.
 * @param[in] M Number of rows of the original matrix.
 * @param[in] N Number of columns of the matrix.
 * @param[in] R The matrix that will hold the result for R. Rank 0 holds the final result.
 * @param[in] rank The rank of the calling process.
 * @param[in] nprocs The number of processes.
 * @param[in] comm The MPI communicator
 */
void TSQR(double *A, const int M, const int N, double *R, const int rank, const int nprocs, MPI_Comm comm){
    
    const int steps = log2(nprocs);
    int start, end;
    decomp1d(M, nprocs, rank, &start, &end); /* Partition M rows over processes */
    // if(!rank){printf("Start: %d, end: %d\n", start, end);}
    
    double * tempA = malloc((end-start+1)*N*sizeof(double)); /* Allocate space to not overwrite A */
    memcpy(tempA, A, (end-start+1)*N*sizeof(double));
    double * tau = malloc(N*sizeof(double)); /* Allocate space to hold tau's */

    for(int step=0;step<=steps;step++){ /* Parallel TSQR loop */

        if(is_active(rank, step)){

            lapack_int rows = ((!step)?end-start+1 : 2*N); 
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
            printf("Rank %d is here\n", rank);
            if(step<steps){
                if(is_active(rank, step + 1)){
                    printf("here!");
                    /* Receive R from other process */
                    if(realloc(tempA, 2*N*N*sizeof(double)) == NULL){
                        perror("Could not reallocate memory in tsqr");
                    }
                    memcpy(tempA, R, N*N*sizeof(double));
                    MPI_Recv(tempA + N*N, N*N, MPI_DOUBLE, MPI_ANY_SOURCE, 1, comm, MPI_STATUS_IGNORE);
                }else{
                    /* Send R to other active process */
                    int lower_active = find_lower_active(rank, step + 1);
                    printf("Lower active for rank %d in step %d: %d\n", rank, step, lower_active);
                    MPI_Send(R, N*N, MPI_DOUBLE, lower_active, 1, comm);
                }
            }
        }else{
            break;
        }
    }
    free(tau);
    free(tempA);
}