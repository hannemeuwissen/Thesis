/**
 * @file sparse.c
 * @brief Code related to working with sparse CSR matrices, part of Thesis 
 * project in High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 2.0
 * @date 2023-05-27
 */
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<mkl_types.h>
#include<mkl_cblas.h>
#include"sparse.h"

/**
 * @brief Function that reads a CSR matrix from a txt file. 
 * @param M The sparse_CSR matrix object that will hold the matrix.
 * @param filename The name of the file. This file contains first the number of columns, the number of
 * rows and the number of nonzeros (each on one line), after this it contains the array of row pointers
 * (each value on one line), and then the column indices and values (one index and corresponding value
 * on one line).
 */
void read_CSR(sparse_CSR * M, const char *const filename){
    FILE *fp;
    if((fp = fopen(filename, "r"))==NULL){
		perror("Error trying to open the file that contains the csr matrix.");
		exit(-1);
    }
    int temp;
    if(fscanf(fp, "%d", &temp) == 0){
        perror("Invalid CSR file");
        exit(-1);
    }
    M->ncols = temp;
    if(fscanf(fp, "%d", &temp) == 0){
        perror("Invalid CSR file");
        exit(-1);
    }
    M->nrows = temp;
    if(fscanf(fp, "%d", &temp) == 0){
        perror("Invalid CSR file");
        exit(-1);
    }
    M->nnz = temp;
    M->rowptrs = malloc((M->nrows + 1)*sizeof(int));
    M->colindex = malloc((M->nnz)*sizeof(int));
    M->values = malloc((M->nnz)*sizeof(double));
    for(int i = 0;i<=M->nrows;i++){
        if(fscanf(fp, "%d", M->rowptrs + i ) == 0){
            perror("Incorrect CSR matrix dimensions in file.");
            exit(-1);
        }
    }
    for (int i=0;i<M->nnz;i++){
        if(fscanf(fp, "%d", M->colindex + i) == 0){
            perror("Incorrect CSR matrix dimensions in file.");
            exit(-1);
        }
        if(fscanf(fp, "%lf", M->values + i) == 0){
            perror("Incorrect CSR matrix dimensions in file.");
            exit(-1);
        }
	}
    fclose(fp);
}

/**
 * @brief Function that prints a sparse CSR matrix structure.
 * @param M The sparse_CSR structure to be printed.
 */
void print_CSR(sparse_CSR * M){
    printf("Number of columns and rows: %d\nNumber of nonzeros: %d\n", M->ncols, M->nnz);
    printf("Row pointers: ");
    for(int i=0;i<(M->nrows+1);i++){
        printf("%d ", M->rowptrs[i]);
    }
    printf("\nColumn indices: ");
    for(int i=0;i<(M->nnz);i++){
        printf("%d ", M->colindex[i]);
    }
    printf("\nValues: ");
    for(int i=0;i<(M->nnz);i++){
        printf("%lf ", M->values[i]);
    }
    printf("\n");
}

/**
 * @brief Function that calculates the sparse matrix mector multiplication between
 * a sparse_CSR matrix and a vector with compatible length.
 * @param M Sparse CSR matrix structure.
 * @param v Vector.
 * @param len Length of the vector.
 * @param result Result vector.
 */
void spmv(sparse_CSR M, double * v, double len, double * result){
    if(len != M.nrows){
        perror("incompatible dimensions in spmv.\n");
        exit(EXIT_FAILURE);
    }
    for(int i=0;i<len;i++){
        result[i] = 0.0;
        for(int j=M.rowptrs[i];j<M.rowptrs[i+1];j++){
            result[i] += M.values[j]*v[M.colindex[j]];
        }
    }
}

/**
 * @brief Function that prints a vector of given length.
 * @param v Vector.
 * @param len Length of the vector.
 */
void print_vector(double * v, const int len){
    for(int i = 0;i<len;i++){
        printf("%lf\n", v[i]);
    }
}

/**
 * @brief Function that prints a matrix of given dimensions.
 * @param A The matrix.
 * @param n The number of rows.
 * @param m The number of columns.
 */
void print_matrix(double * A, const int n, const int m){
    for(int i = 0;i<n;i++){
        for(int j=0;j<m;j++){
            printf("%lf ", A[i*m + j]);
        }
        printf("\n");
    }
}

/**
 * @brief Function that reads data from a file and stores it in a matrix.
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

/**
 * @brief Function that computes the (m+1)-Krylov subspace of the sparse CSR matrix A,
 * thus the space spanned by the vectors {b, Ab, ..., A^m b}.
 * @param A Sparse sparse_CSR matrix struct. A should represent a len x len matrix.
 * @param b Vector b.
 * @param len Length of the vector.
 * @param Q Matrix which will hold the resulting basis vectors of the Krylov subspace. (len x m+1)
 * @param H Matrix which will hold the resulting upper Hessenberg matrix H. ((m+1) x m)
 * @param m The number of iterations.
 */
void Arnoldi(sparse_CSR A, double * b, const int len, double * Q, double * H, const int m){
    if(A.ncols != len){
        perror("Incompatible dimensions in Arnoldi.");
        exit(EXIT_FAILURE);
    }
    double eps = 1e-12; /* Tolerance for h: stop if h too small */
    double * Q_trans = malloc(len*(m+1)*sizeof(double));
    cblas_dcopy(len, b, 1, Q_trans, 1); /* Set first vector of Q */
    double norm_b = cblas_dnrm2(len, Q_trans, 1);
    cblas_dscal(len, 1/norm_b, Q_trans, 1); /* Normalize */
    double * w = malloc(len*sizeof(double));
    for(int j=1;j < m+1;j++){
        spmv(A, Q_trans + (j-1)*len, len, w);
        for(int i=0;i<j;i++){
            H[i*m + (j-1)] = cblas_ddot(len, w, 1, Q_trans + i*len, 1);
            cblas_daxpy(len,-H[i*m + (j-1)],Q_trans + i*len,1,w, 1);
        }
        H[j*m + (j-1)] = cblas_dnrm2(len, w, 1);
        if(H[j*m + (j-1)] < eps){
            break;
        }
        cblas_dscal(len, 1/H[j*m + (j-1)], w, 1);
        cblas_dcopy(len, w, 1, Q_trans + j*len, 1);
    }
    /* Transpose final result */
    for(int i=0;i<len;i++){
        for(int j=0;j<m+1;j++){
            Q[i*(m+1) + j] = Q_trans[j*len + i];
        }
    }
}