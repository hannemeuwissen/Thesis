/**
 * @file tsqr.c
 * @brief File containing function definitions related to communication avoiding
 * tall skinny QR decomposition for case 1 for MAP55672 Case Studies for HPC at Trinity College Dublin.
 * @author Hanne Meuwissen (22307813)
 * @version 2.0
 * @date 2023-02-13
 */

#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<string.h>
#include<mkl_types.h>
#include<mkl_cblas.h>
#include<math.h>
#include"mkl.h"
#include"matrix.h"

/**
 * @brief Function that performs the communication avoiding TSQR using OpenMP.
 * @param[in] A The matrix that will hold the data. It will hold the final result for Q.
 * @param[in] M Number of rows of the matrix.
 * @param[in] N Number of columns of the matrix.
 * @param[in] R The matrix that will hold the result for R.
 * @param[in] n The number of processes to use.
 */
void TSQR(double *A, const int M, const int N, double *R, const int n){
    const int steps = log2(n);
    
    double * block, * tau, *tempA;
    int first_partition_s[n];
    int first_partition_e[n];
    if(n>1){ /* Allocate memory for partial product to form final Q */
        block = malloc(M*2*N*sizeof(double));
    }

    for(int step=0;step<=steps;step++){ /* Parallel TSQR loop */

        int active_procs = n/pow(2.0, (double) step); /* The number of processes active in this step */
        tau = malloc(N*active_procs*sizeof(double)); /* Allocate space to hold tau's */
        tempA = malloc(N*N*active_procs*sizeof(double)); /* Allocate space to temporarily hold local R's */

        #pragma omp parallel num_threads(active_procs)
        {
            int rank = omp_get_thread_num();
            int start,end;
            if(!step){ 
                decomp1d(M, active_procs, rank, &start, &end); /* Partition M rows over processes */
                first_partition_s[rank] = start; /* Save the partition */
                first_partition_e[rank] = end;
            }else{
                start = rank*2*N;
                end = start+2*N-1;
            }
            lapack_int rows = end-start+1; 
            lapack_int cols = N;

            /* Calculate local Housholder QR in terms of tau and v's */
            LAPACKE_dgeqrf(CblasRowMajor, rows, cols, A + start*N, cols, tau + rank*N); 
            
            /* Save R parts in temporary space for next iteration */
            for(int i=0;i<N;i++){ 
                for(int j=0;j<N;j++){
                    tempA[j + (rank*N + i)*N] = ((i>j) ? 0 : A[j + (start + i)*N]);
                }
            }
            
            /* Get local Q */
            LAPACKE_dorgqr(CblasRowMajor, rows, cols, cols, A + start*N, cols, tau + rank*N); 
            
            if(step>0){ /* Perform part of multiplication for final Q */
                int skip = n/active_procs;
                int qstart = first_partition_s[rank*skip]; /* Determine starting point based on partition of M rows */
                int qend = first_partition_e[(rank+1)*skip -1];
                int qrows = qend-qstart+1;
                double * res = malloc(qrows*N*sizeof(double));
                /* Perform local part of product */
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, qrows, N, 2*N, 1.0, block + qstart*2*N, 2*N, A + start*N, N, 0.0, res, N);             
                if(step<steps){ /* Place in next block */
                    memset(block + qstart*2*N, 0, qrows*2*N*sizeof(double));
                    for(int i=0;i<qrows;i++){
                        for(int j=0;j<N;j++){
                            block[j + N*(rank%2) + (qstart + i)*2*N] = res[j + i*N];
                        }
                    }
                }else{ /* Replace A with last result */
                    set_equal(A, res, M, N);
                }
                free(res);
            }else if((step == 0) && (steps!=0)){ /* Place in block */
                memset(block + start*2*N, 0, rows*2*N*sizeof(double));
                for(int i=0;i<rows;i++){
                    for(int j=0;j<N;j++){
                        block[j + N*(rank%2) + (start + i)*2*N] = A[j + (start + i)*N];
                    }
                }
            }
        }
        free(tau);
        
        if(step == steps){ /* Get final R from temporary space */
            for(int i=0;i<N;i++){
                for(int j=0;j<N;j++){
                    R[j + i*N] = ((i>j) ? 0 : tempA[j + i*N]);
                }
            }
        }else{ /* Place R's in A for next local QR step */
            set_equal(A, tempA, N*active_procs, N);
            free(tempA);
        }
    }
    if(steps>0){
        free(block);
    }
}