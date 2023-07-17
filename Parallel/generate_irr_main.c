/**
 * @file main.c
 * @brief Main code to investigate CA-Arnoldi, part of Thesis project in 
 * High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 4.0
 * @date 2023-06-02
 */
#include<stdlib.h>
#include<stdio.h>
#include"graph.h"


/**
 * @brief Function that saves a CSR matrix to a file.
 * @param filename_A The name of the file.
 * @param A The CSR_matrix.
 */
void save_CSR(char * filename_A, sparse_CSR A){
    FILE *fp;
    fp = fopen(filename_A, "w");
    if(!fp){
      fprintf(stderr, "Error: can't open file %s\n",filename_A);
      exit(4);
    }
    fprintf(fp, "%d\n", A.ncols);
    fprintf(fp, "%d\n", A.nrows);
    fprintf(fp, "%d\n", A.nnz);
    for(int i=0;i<(A.nrows+1);i++){
        fprintf(fp, "%d\n", A.rowptrs[i]);
    }
    for(int i=0;i<A.nnz;i++){
        fprintf(fp,"%d %lf\n", A.colindex[i], A.values[i]);
    }
    fclose(fp);
}

int main(int argc, char **argv){

    if(argc < 5){
        fprintf(stderr, "Usage: %s file M min_nnz max_nnz\n", argv[0]);
        return 1;
    }

    int M, min_nnz, max_nnz;
    if((sscanf(argv[2],"%d",&M) == 0) || (sscanf(argv[3],"%d",&min_nnz) == 0) || (sscanf(argv[4],"%d",&max_nnz) == 0)){
        fprintf(stderr, "Usage: %s file M min_nnz max_nnz\n", argv[0]);
        return 1;
    }

    sparse_CSR A = generate_irregular_csr(M, min_nnz, max_nnz);
    save_CSR(argv[1], A);
    
    return 0;
}