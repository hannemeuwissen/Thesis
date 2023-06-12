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

/**
 * @brief Function that reads and parses command line arguments.
 * @param[in] argc The command line argument count.
 * @param[in] argv Pointer to the command line arguments.
 * @param[in] n Integer pointer to nr of rows.
 * @param[in] m Integer pointer to nr of cols.
 * @param[in] seed Integer pointer to seed.
 * @param[in] time Integer pointer to time specifier.
 */
void parse_command_line(const int argc, char * const *argv, char * filename_A, int * M, int * N, int * degree, int * s){
    int c=0;
    while((c = getopt(argc, argv, "d:s:f:m:n:")) != -1){
        switch(c){
            case 'f':
                if(sscanf(optarg,"%s",filename) == 0){
                    printf("Usage : ./main [-f filename] [-m nr of rows] [-n nr of columns] [-d degree] [-s blocksize]\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                break;
            case 'm':
                if(sscanf(optarg,"%d",M) == 0){
                    printf("Usage : ./main [-f filename] [-m nr of rows] [-n nr of columns] [-d degree] [-s blocksize]\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                break;
            case 'n':
                if(sscanf(optarg,"%d",N) == 0){
                    printf("Usage : ./main [-f filename] [-m nr of rows] [-n nr of columns] [-d degree] [-s blocksize]\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                break;
            case 'd':
                if(sscanf(optarg,"%d",degree) == 0){
                    printf("Usage : ./main [-f filename] [-m nr of rows] [-n nr of columns] [-d degree] [-s blocksize]\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                break;
            case 's':
                if(sscanf(optarg,"%d",s) == 0){
                    printf("Usage : ./main [-f filename] [-m nr of rows] [-n nr of columns] [-d degree] [-s blocksize]\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                break;
            case '?':
                printf("Usage : ./main [-f filename] [-m nr of rows] [-n nr of columns] [-d degree] [-s blocksize]\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    if(degree%s != 0){
        printf("Invalid input: the degree of the Krylov subspace should be a multiple of the blocksize (s)\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

int main(int argc, char **argv)
{  
    int myid, nprocs;
    int degree,s,M,N;
    char filename_A[100];

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
    if(!myid){
        parse_command_line(argc, argv, filename_A, &M, &N, &degree, &s);
    }
    MPI_Bcast(&degree, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(filename_A, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int steps = degree/s;
    int start, end;
    decomp1d(M, nprocs, myid, &start, &end);
    int m = end - start + 1;

    // For all blocks:
    for(int block = 0;block < steps;block++){
        // 1. Matrix powers kernel using parallel spmv for size of block
        // NOTE: since output for each spmv is vector --> need to transpose before next part!

        // 2. Block-GS to orthogonalize compared to previous blocks (not the first time)
        bgs(V, W, m, block*s + 1, s, MPI_COMM_WORLD);

        // 3. Orthogonalize block using parallel CA-TSQR
        double * R = malloc(s*s*sizeof(double)); 
        TSQR(A, m, s, R, myid, nprocs, MPI_COMM_WORLD); // Is R needed?
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    return 0;
}