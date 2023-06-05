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

int main(int argc, char **argv)
{  
    int myid, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Read input: degree of Krylov subspace

    // Decide on s in s-step process

    // For all blocks:
    // 1. Matrix powers kernel using parallel spmv 
    // 2. Block-GS to orthogonalize compared to previous blocks (not the first time)
    // 3. Orthogonalize block using parallel CA-TSQR
    
    MPI_Finalize();
    return 0;
}