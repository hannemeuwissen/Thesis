/**
 * @file angle_subspaces.c
 * @brief Main code to investigate the stability of CA-Arnoldi with respect to the block
 * size, scipt calculates the angle between two subspaces; part of Thesis project in High 
 * Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 1.0
 * @date 2023-07-24
 */
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<mkl_cblas.h>
#include<mkl_types.h>
#include"matrix.h"
#include"mpi.h"

/**
 * @brief Function that calculates the Grassmann distance based on singular values.
 * @param s The singular values.
 * @param n The number of singular values.
 * @return The Grassmann distance.
 */
double grassmann_distance(double *s, const int n){
    double res = 0;
    double term;
    for(int i=0;i<n;i++){
        term = acos(s[i]);
        res += pow(term,2);
    }
    return sqrt(res);
}

int main(int argc, char **argv){
    
    if(argc < 5){
        fprintf(stderr, "Usage: %s fileQ1 fileQ2 M n\n", argv[0]);
        return 1;
    }

    int M, n;
    if((sscanf(argv[3],"%d",&M) == 0) || (sscanf(argv[4],"%d",&n) == 0)){
        fprintf(stderr, "Usage: %s fileQ1 fileQ2 M n\n", argv[0]);
        return 1;
    }

    /* Read Q1, Q2 */
    double * Q1 = malloc(M*n*sizeof(double));
    read_matrix_from_file_double(argv[1], 0, Q1, M, n);
    double * Q2 = malloc(M*n*sizeof(double));
    read_matrix_from_file_double(argv[2], 0, Q2, M, n);

    /* Calculate Q1^TQ2 */
    double * D = malloc(n*n*sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, M, 1.0, Q1, n, Q2, n, 0.0, D, n);
    
    /* Singular value decomposition of D */
    double * S = malloc(n*sizeof(double));
    double U,V;
    double * superb = malloc((n-2)*sizeof(double));
    int ret = LAPACKE_dgesvd(CblasRowMajor, "N", "N", n, n, D, n, S, U, 1, V, 1, superb);
    if(ret!=0){
        if(ret<0){
            fprintf(stderr, "LAPACKE_dgesvd failed. Parameter %d had an illegal value\n", abs(ret));
            exit(EXIT_FAILURE);
        }
        else{
            fprintf(stderr, "LAPACKE_dgesvd failed. Leading minor of order %d is not positive definite\n", ret);
            exit(EXIT_FAILURE);
        }
    }

    /* Calculate distance */
    double dist = grassmann_distance(S, n);

    printf("The Grassmann distance between the subspaces in %s and %s is:\n%lf\n", argv[1], argv[2], dist);

    return 0;
}