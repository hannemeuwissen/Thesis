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
#include"hess.h"
#include"mkl.h"

/**
 * @brief Function that handles the input arguments. 
 * @param argc
 * @param argv 
 * @param filename_A The name of the file that contains the sparse CSR matrix.
 * @param M The number of rows.
 * @param N The number of columns.
 * @param filename_v The name of the file that contains the initial vector.
 * @param degree The degree of the Krylov subspace.
 * @param s The blocksize (s-step Krylov subspace method).
 */
void parse_command_line(const int argc, char * const *argv, char * filename_A, int * M, int * N, char * filename_v, int * degree, int * s){
    int c=0;
    while((c = getopt(argc, argv, "d:s:f:v:m:n:")) != -1){
        switch(c){
            case 'f':
                if(sscanf(optarg,"%s",filename_A) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of rows] [-n nr of columns] [-d degree] [-s blocksize]\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                break;
            case 'v':
                if(sscanf(optarg,"%s",filename_v) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of rows] [-n nr of columns] [-d degree] [-s blocksize]\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                break;
            case 'm':
                if(sscanf(optarg,"%d",M) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of rows] [-n nr of columns] [-d degree] [-s blocksize]\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                break;
            case 'n':
                if(sscanf(optarg,"%d",N) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of rows] [-n nr of columns] [-d degree] [-s blocksize]\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                break;
            case 'd':
                if(sscanf(optarg,"%d",degree) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of rows] [-n nr of columns] [-d degree] [-s blocksize]\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                break;
            case 's':
                if(sscanf(optarg,"%d",s) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of rows] [-n nr of columns] [-d degree] [-s blocksize]\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                break;
            case '?':
                printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of rows] [-n nr of columns] [-d degree] [-s blocksize]\n");
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
    char filename_A[100], filename_v[100];

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

    // Read input: degree of Krylov subspace + s (-> degree multiple of s) + A + size A (M,N) + input vector
    if(!myid){
        parse_command_line(argc, argv, filename_A, filename_v, &M, &N, &degree, &s);
    }
    MPI_Bcast(&degree, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(filename_A, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(filename_v, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
    int steps = degree/s;
    int start, end;
    decomp1d(M, nprocs, myid, &start, &end);
    int m = end - start + 1;

    // Read data

    // For all blocks:
    for(int block = 0;block < steps;block++){

        // Fix B

        // Matrix powers kernel using parallel spmv for size of block
        // NOTE: since output for each spmv is vector --> need to transpose before next part!

        if(!block){
            // Orthogonalize block using parallel CA-TSQR
            double * R_ = malloc((s+1)*(s+1)*sizeof(double)); 
            TSQR(V, m, s, R_, myid, nprocs, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            // set mathcal Q
            double * mathcalQ = malloc(m*(s+1)*sizeof(double));
            memcpy(mathcalQ, V, m*(s+1)*sizeof(double));

            // Calculate mathcal H
            if(!myid){
                double * mathcalH = malloc((s+1)*s*sizeof(double));

                // set B_

                double * R = malloc(s*s*sizeof(double)); 
                get_R(R, R_, s+1);
                calc_hess(mathcalH, R_, B_, R, s+1, s);
            }
        }else{
            // Block-CGS to orthogonalize compared to previous blocks
            double mathcalR = malloc((block*s + 1)*s*sizeof(double));
            bgs(mathcalQ, V, mathcalR, m, block*s + 1, s, MPI_COMM_WORLD);

            // Orthogonalize block using parallel CA-TSQR
            double * R = malloc(s*s*sizeof(double)); 
            TSQR(V, m, s, R, myid, nprocs, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD); 

            // Update mathcal Q

            // Update mathcal R_k

            // Update mathcal B

            // Update mathcal H
        }

    }
    
    MPI_Finalize();
    return 0;
}