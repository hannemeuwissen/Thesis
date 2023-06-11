/**
 * @file main.c
 * @author Hanne Meuwissen (22307813)
 * @brief Code for thesis at Trinity College Dublin.
 * @version 1.0
 * @date 2023-06-02
 */
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include<mkl_cblas.h>
#include<mkl_types.h>
#include"matrix.h"
#include"graph.h"
#include"sparse.h"
#include"tsqr_mpi.h"
#include"bgs.h"
#include"mkl.h"

int main(int argc, char **argv)
{  
    int myid, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(!myid){
        float logprocs = log2(nprocs);
        if(ceil(logprocs) != floor(logprocs)){
            fprintf(stderr,"Error: The number of processes needs to be a power of 2 (because of TSQR).\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Read input: degree of Krylov subspace + s (-> degree multiple of s) + A + size A (M,N)
    int steps = degree/s;
    int start, end;
    decomp1d(M, nprocs, myid, &start, &end);
    int m = end - start + 1;

    // For all blocks:
    
    // 1. Matrix powers kernel using parallel spmv for size of block
    // NOTE: since output for each spmv is vector --> need to transpose in next part!

    // 2. Block-GS to orthogonalize compared to previous blocks (not the first time)
    // First transpose!
    bgs(V, W, m, N, MPI_COMM_WORLD);

    // 3. Orthogonalize block using parallel CA-TSQR
    double * R = malloc(n*n*sizeof(double)); 
    TSQR(A, m*nprocs, n, R, myid, nprocs, MPI_COMM_WORLD); // Is R needed?
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Finalize();
    return 0;
}