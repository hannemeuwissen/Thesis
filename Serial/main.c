/**
 * @file main.c
 * @brief Main code for testing serial Arnoldi, part of Thesis 
 * project in High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 1.0
 * @date 2023-05-26
 */
#include<stdlib.h>                                                             
#include<stdio.h>
#include<string.h>
#include"graph.h"
#include"sparse.h"                                                       
                                                                                
int main(int argc, char **argv){
    int n,m;

    if(argc < 4){
        fprintf(stderr, "Usage: %s filename_A filename_v degree\n", argv[0]);
        return 1;
    }
    if((sscanf(argv[3],"%d",&m) == 0)){
        fprintf(stderr, "Usage: %s filename_A filename_v degree\n", argv[0]);
        return 1;
    }

    /* Check serial arnoldi */

    sparse_CSR T;
    read_CSR(&T, argv[1]);
    int n = T.nrows;

    double * Q = malloc(n*m*sizeof(double));
    double * H = malloc(m*(m-1)*sizeof(double));
    memset(H, 0, m*(m-1)*sizeof(double));
    double *b = malloc(n*sizeof(double));
    read_matrix_from_file(argv[2], 0, b, n, 1);
    
    Arnoldi(T, b, n, Q, H, m);
    print_matrix(Q, n, m);
    print_matrix(H, m, m-1);
    free(H);
    free(Q);
    free(b);
}     