/**
 * @file parse.c
 * @author Hanne Meuwissen (22307813)
 * @brief Code for thesis at Trinity College Dublin.
 * @version 1.0
 * @date 2023-06-02
 */
#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include<mpi.h>

/**
 * @brief Function that handles the input arguments. 
 * @param argc
 * @param argv
 * @param M The number of nodes.
 * @param nnz The number of edges per node in the random generated graph.
 * @param filename_v The name of the file that contains the initial vector.
 * @param degree The degree of the Krylov subspace.
 * @param s The blocksize (s-step Krylov subspace method).
 */
void parse_command_line_regular(const int argc, char * const *argv, int * M, int * nnz, char * filename_v, int * degree, int * s, MPI_Comm comm){
    int c=0;
    *M = 800;
    *nnz = 50;
    while((c = getopt(argc, argv, "d:s:v:m:n:z:")) != -1){
        switch(c){
            case 'v':
                if(sscanf(optarg,"%s",filename_v) == 0){
                    printf("Usage : ./main [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 'm':
                if(sscanf(optarg,"%d",M) == 0){
                    printf("Usage : ./main [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 'n':
                if(sscanf(optarg,"%d",N) == 0){
                    printf("Usage : ./main [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 'z':
                if(sscanf(optarg,"%d",nnz) == 0){
                    printf("Usage : ./main [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 'd':
                if(sscanf(optarg,"%d",degree) == 0){
                    printf("Usage : ./main [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 's':
                if(sscanf(optarg,"%d",s) == 0){
                    printf("Usage : ./main [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case '?':
                printf("Usage : ./main [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize]\n");
                MPI_Abort(comm, 1);
        }
    }
}

/**
 * @brief Function that handles the input arguments. 
 * @param argc
 * @param argv 
 * @param filename_A The name of the file that contains the sparse CSR matrix.
 * @param M The number of nodes.
 * @param filename_v The name of the file that contains the initial vector.
 * @param degree The degree of the Krylov subspace.
 * @param s The blocksize (s-step Krylov subspace method).
 */
void parse_command_line_irregular(const int argc, char * const *argv, char * filename_A, int * M, char * filename_v, int * degree, int * s, MPI_Comm comm){
    int c=0;
    *M = 800;
    while((c = getopt(argc, argv, "d:s:f:v:m:n:")) != -1){
        switch(c){
            case 'f':
                if(sscanf(optarg,"%s",filename_A) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 'v':
                if(sscanf(optarg,"%s",filename_v) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 'm':
                if(sscanf(optarg,"%d",M) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 'n':
                if(sscanf(optarg,"%d",N) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 'd':
                if(sscanf(optarg,"%d",degree) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 's':
                if(sscanf(optarg,"%d",s) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case '?':
                printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-d degree] [-s blocksize]\n");
                MPI_Abort(comm, 1);
        }
    }
}