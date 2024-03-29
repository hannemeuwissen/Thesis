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
#include<string.h>
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
 * @param t Indicator to print out the timings.
 * @param comm The MPI communicator.
 */
void parse_command_line_regular(const int argc, char * const *argv, int * M, int * nnz, char * filename_v, int * degree, int * s, int * t, MPI_Comm comm){
    int c=0;
    *M = 800;
    *nnz = 50;
    strncpy(filename_v, "v.txt", 100);
    *degree = 8;
    *s = 2;
    *t = 0;
    while((c = getopt(argc, argv, "ad:s:v:m:z:t")) != -1){
        switch(c){
            case 'a':
                printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                MPI_Abort(comm, 1);
                break;
            case 'v':
                if(sscanf(optarg,"%s",filename_v) == 0){
                    printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 'm':
                if(sscanf(optarg,"%d",M) == 0){
                    printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 'z':
                if(sscanf(optarg,"%d",nnz) == 0){
                    printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 'd':
                if(sscanf(optarg,"%d",degree) == 0){
                    printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 's':
                if(sscanf(optarg,"%d",s) == 0){
                    printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 't':
                *t = 1;
                break;
            case '?':
                printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                MPI_Abort(comm, 1);
        }
    }
}

/**
 * @brief Function that handles the input arguments.
 * @param argc First standard input parameter.
 * @param argv Second standard input parameter.
 * @param M The number of nodes in the full graph.
 * @param min_nnz The minimum number of edges per node in the random generated graph.
 * @param max_nnz The maximum number of edges per node in the random generated graph.
 * @param filename_v The name of the file that contains the initial vector.
 * @param degree The degree of the Krylov subspace.
 * @param s The blocksize (s-step Krylov subspace method).
 * @param t Indicator to print out the timings.
 * @param comm The MPI communicator.
 */
void parse_command_line_irregular(const int argc, char * const *argv, int * M, int * min_nnz, int * max_nnz, char * filename_v, int * degree, int * s, int * t, MPI_Comm comm){
    int c=0;
    *M = 800;
    *min_nnz = 0;
    *max_nnz = 800;
    strncpy(filename_v, "v.txt", 100);
    *degree = 8;
    *s = 2;
    *t = 0;
    while((c = getopt(argc, argv, "ad:s:tv:m:z:x:")) != -1){
        switch(c){
            case 'a':
                printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                MPI_Abort(comm, 1);
                break;
            case 'v':
                if(sscanf(optarg,"%s",filename_v) == 0){
                    printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 'z':
                if(sscanf(optarg,"%d",min_nnz) == 0){
                    printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 'x':
                if(sscanf(optarg,"%d",max_nnz) == 0){
                    printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 'm':
                if(sscanf(optarg,"%d",M) == 0){
                    printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 'd':
                if(sscanf(optarg,"%d",degree) == 0){
                    printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 's':
                if(sscanf(optarg,"%d",s) == 0){
                    printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 't':
                *t = 1;
                break;
            case '?':
                printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z min nnz] [-x max nnz] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                MPI_Abort(comm, 1);
        }
    }
}


/**
 * @brief Function that handles the input arguments.
 * @param argc First standard input parameter.
 * @param argv Second standard input parameter.
 * @param filename_A The name of the file that contains the sparse CSR matrix.
 * @param filename_v The name of the file that contains the initial vector.
 * @param degree The degree of the Krylov subspace.
 * @param s The blocksize (s-step Krylov subspace method).
 * @param t Indicator to print out the timings.
 * @param lb Indicator if load-balancing should be added or not.
 * @param q Indicator if Q needs to be saved to a file.
 * @param h Indicator if H needs to be saved to a file.
 * @param comm The MPI communicator.
 */
void parse_command_line_lb(const int argc, char * const *argv, char * filename_A, char * filename_v, int * degree, int * s, int * t, int * lb, int * q, int * h, MPI_Comm comm){
    int c=0;
    strncpy(filename_A, "test100irr.txt", 100);
    strncpy(filename_v, "v.txt", 100);
    *lb = 0;
    *degree = 8;
    *s = 2;
    *q = 0;
    *h = 0;
    *t = 0;
    while((c = getopt(argc, argv, "ad:f:s:tv:lqh")) != -1){
        switch(c){
            case 'a':
                printf("Usage : %s [-v filename_v] [-m nr of nodes] [-z nnz per row] [-d degree] [-s blocksize] [-t] [-a]\n", argv[0]);
                MPI_Abort(comm, 1);
                break;
            case 'f':
                if(sscanf(optarg,"%s",filename_A) == 0){
                    printf("Usage : %s [-f filename_A] [-v filename_v] [-d degree] [-s blocksize] [-t] [-a] [-l] [-q] [-h]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 'v':
                if(sscanf(optarg,"%s",filename_v) == 0){
                    printf("Usage : %s [-f filename_A] [-v filename_v] [-d degree] [-s blocksize] [-t] [-a] [-l] [-q] [-h]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 'd':
                if(sscanf(optarg,"%d",degree) == 0){
                    printf("Usage : %s [-f filename_A] [-v filename_v] [-d degree] [-s blocksize] [-t] [-a] [-l] [-q] [-h]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 's':
                if(sscanf(optarg,"%d",s) == 0){
                    printf("Usage : %s [-f filename_A] [-v filename_v] [-d degree] [-s blocksize] [-t] [-a] [-l] [-q] [-h]\n", argv[0]);
                    MPI_Abort(comm, 1);
                }
                break;
            case 't':
                *t = 1;
                break;
            case 'l':
                *lb = 1;
                break;
            case 'q':
                *q = 1;
                break;
            case 'h':
                *h = 1;
                break;
            case '?':
                printf("Usage : %s [-f filename_A] [-v filename_v] [-d degree] [-s blocksize] [-t] [-a] [-l] [-q] [-h]\n", argv[0]);
                MPI_Abort(comm, 1);
        }
    }
}