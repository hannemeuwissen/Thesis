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
#include"mkl.h"
#include<mkl_cblas.h>
#include<mkl_types.h>


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
        term = ((s[i] > 1) ? 1 : s[i]);
        res += pow(acos(term),2);
    }
    return sqrt(res);
}

/**
 * @brief Function that reads data from a file (double precision) and stores it in a matrix.
 * @param[in] filename File to read the matrix from.
 * @param[in] skip The number of elements to skip.
 * @param[in] A The matrix that will hold the data.
 * @param[in] M Number of rows of the matrix.
 * @param[in] N Number of columns of the matrix.
 */
void read_matrix_from_file(const char *const filename, const int skip, double *A, const int M, const int N){
    FILE *fp;
    if((fp = fopen(filename, "r"))==NULL){
		perror("Error trying to open the file");
		exit(-1);
    }
    double temp;
    for(int r=0;r<skip;r++){
        if(fscanf(fp, "%lf", &temp) == 0){
            perror("Incorrect matrix dimensions");
            exit(-1);
        }
    }
    for (int i=0;i<M;i++){
        for(int j=0;j<N;j++){
            if(fscanf(fp, "%lf", A + j + i*N) == 0){
                printf("Goes wrong on row %d and col %d\n", i, j);
			    perror("Incorrect matrix dimensions");
			    exit(-1);
		    }
        }
	}
    fclose(fp);
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
    read_matrix_from_file(argv[1], 0, Q1, M, n);
    double * Q2 = malloc(M*n*sizeof(double));
    read_matrix_from_file(argv[2], 0, Q2, M, n);

    /* Calculate Q1^TQ2 */
    double * D = malloc(n*n*sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, M, 1.0, Q1, n, Q2, n, 0.0, D, n);
    
    /* Singular value decomposition of D */
    double * S = malloc(n*sizeof(double));
    double U,V;
    double * superb;
    int ret = LAPACKE_dgesvd(CblasRowMajor, 'N', 'N', n, n, D, n, S, &U, 1, &V, 1, superb);
    if(ret!=0){
        if(ret<0){
            fprintf(stderr, "LAPACKE_dgesvd failed. Parameter %d had an illegal value\n", abs(ret));
            exit(EXIT_FAILURE);
        }
    }

    /* Calculate distance */
    double dist = grassmann_distance(S, n);

    printf("The Grassmann distance between the subspaces in %s and %s is:\n%e\n", argv[1], argv[2], dist);

    free(Q1);
    free(Q2);
    free(S);
    free(D);

    return 0;
}