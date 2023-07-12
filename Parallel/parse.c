/**
 * @file parse.c
 * @brief Code related to parsing the input parameters to investigate CA-Arnoldi, part of
 * Thesis project in High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 2.0
 * @date 2023-06-02
 */
#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include<mpi.h>

/**
 * @brief Function that handles the input arguments. 
 * @param argc First standard input parameter.
 * @param argv Second standard input parameter.
 * @param M The number of nodes in the full graph.
 * @param nnz The number of edges per node in the random generated graph.
 * @param filename_v The name of the file that contains the initial vector.
 * @param degree The degree of the Krylov subspace.
 * @param s The blocksize (s-step Krylov subspace method).
 */
void parse_command_line_regular(const int argc, char * const *argv, int * M, int * nnz, char * filename_v, int * degree, int * s, MPI_Comm comm){
    int c=0;
    *M = 800;
    *nnz = 50;
    while((c = getopt(argc, argv, "d:s:v:m:z:")) != -1){
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
                (*degree)--;
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
 * @param argc First standard input parameter.
 * @param argv Second standard input parameter.
 * @param filename_A The name of the file that contains the sparse CSR matrix.
 * @param M The number of nodes in the full graph.
 * @param min_nnz The minimum number of edges per node in the random generated graph.
 * @param max_nnz The maximum number of edges per node in the random generated graph.
 * @param filename_v The name of the file that contains the initial vector.
 * @param degree The degree of the Krylov subspace.
 * @param s The blocksize (s-step Krylov subspace method).
 */
void parse_command_line_irregular(const int argc, char * const *argv, char * filename_A, int * M, int * min_nnz, int * max_nnz, char * filename_v, int * degree, int * s, MPI_Comm comm){
    int c=0;
    *M = 800;
    *min_nnz = 20;
    *max_nnz = 40;
    while((c = getopt(argc, argv, "d:s:f:v:m:z:x:")) != -1){
        switch(c){
            case 'f':
                if(sscanf(optarg,"%s",filename_A) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 'v':
                if(sscanf(optarg,"%s",filename_v) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 'z':
                if(sscanf(optarg,"%d",min_nnz) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 'x':
                if(sscanf(optarg,"%d",max_nnz) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 'm':
                if(sscanf(optarg,"%d",M) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case 'd':
                if(sscanf(optarg,"%d",degree) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                (*degree)--;
                break;
            case 's':
                if(sscanf(optarg,"%d",s) == 0){
                    printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize]\n");
                    MPI_Abort(comm, 1);
                }
                break;
            case '?':
                printf("Usage : ./main [-f filename_A] [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize]\n");
                MPI_Abort(comm, 1);
        }
    }
}