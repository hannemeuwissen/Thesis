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
#include"graph.h"
#include"sparse.h"
#include"tsqr_mpi.h"

int main(int argc, char **argv)
{  
    int myid, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    // /* Test SPMV */
    // sparse_CSR M = generate_regular_graph_part_csr(2, 8, 3);
    // /* Print in order */
    // if(!myid){
    //     printf("Rank %d:\n", myid);
    //     print_CSR(&M);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(myid == 1){
    //     printf("Rank %d:\n", myid);
    //     print_CSR(&M);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(myid == 2){
    //     printf("Rank %d:\n", myid);
    //     print_CSR(&M);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(myid == 3){
    //     printf("Rank %d:\n", myid);
    //     print_CSR(&M);
    // }
    // double x[2] = {2.0,1.0};
    // double result[2];
    // spmv(M, x, 2, result, myid, nprocs, MPI_COMM_WORLD);
    // if(!myid){
    //     print_vector(result, 2); // result should be 1 overall (sum of row elements)
    // }

    /* Test TSQR */
    // let each process read the same matrix part
    double A[20] = {9, 2, 7, 5, 10, 9, 8, 10, 9, 9, 8, 3, 3, 0, 4, 10, 6, 7, 10, 0};
    double * R = malloc(4*4*sizeof(double));
    TSQR(A, 20, 4, R, myid, nprocs, MPI_COMM_WORLD);
    if(!myid){
        print_matrix(R, 4, 4);
    }
    // check resulting R
    // do lapack algorithm to get Q and check this

    // Read input: degree of Krylov subspace

    // Decide on s in s-step process

    // For all blocks:
    // 1. Matrix powers kernel using parallel spmv 
    // 2. Block-GS to orthogonalize compared to previous blocks (not the first time)
    // 3. Orthogonalize block using parallel CA-TSQR
    
    MPI_Finalize();
    return 0;
}