/**
 * @file matrix.c
 * @brief Code related to working with dense matrices, part of Thesis project in 
 * High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 3.0
 * @date 2023-06-02
 */
#include<stdlib.h>
#include<stdio.h>
#include<mpi.h>

/**
 * @brief Function that calculates the average of a given vector.
 * @param v The vector.
 * @param n The length of the vector.
 * @return The average.
 */
double average(double * v, const int n){
    double res = 0;
    for(int i=0;i<n;i++){
        res += v[i];
    }
    return res/n;
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
 * @brief Function that saves a matrix of given dimensions to a file.
 * @param A The matrix.
 * @param n The number of rows.
 * @param m The number of columns.
 * @param filename The name of the file.
 */
void save_matrix(double * A, const int n, const int m, char * filename){
    FILE *fp;
    fp = fopen(filename, "w");
    if(!fp){
      fprintf(stderr, "Error: can't open file %s\n",filename);
      exit(4);
    }   
    for(int i = 0;i<n;i++){
        for(int j=0;j<m;j++){
           fprintf(fp, "%lf ", A[i*m + j]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

/**
 * @brief Function that saves a matrix of given dimensions to a file (in double precision).
 * @param A The matrix.
 * @param n The number of rows.
 * @param m The number of columns.
 * @param filename The name of the file.
 */
void save_matrix_double(double * A, const int n, const int m, char * filename){
    FILE *fp;
    fp = fopen(filename, "w");
    if(!fp){
      fprintf(stderr, "Error: can't open file %s\n",filename);
      exit(4);
    }   
    for(int i = 0;i<n;i++){
        for(int j=0;j<m;j++){
           fprintf(fp, "%.17g ", A[i*m + j]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

/**
 * @brief Function that prints the transpose of a matrix of given dimensions.
 * @param A The matrix.
 * @param n The number of rows (before transpose).
 * @param m The number of columns (before transpose).
 */
void print_matrix_transposed(double * A, const int n, const int m){
    for(int j = 0;j<m;j++){
        for(int i=0;i<n;i++){
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
 * @brief Function that reads data from a file (double precision) and stores it in a matrix.
 * @param[in] filename File to read the matrix from.
 * @param[in] skip The number of elements to skip.
 * @param[in] A The matrix that will hold the data.
 * @param[in] M Number of rows of the matrix.
 * @param[in] N Number of columns of the matrix.
 */
void read_matrix_from_file_double(const char *const filename, const int skip, double *A, const int M, const int N){
    FILE *fp;
    if((fp = fopen(filename, "r"))==NULL){
		perror("Error trying to open the file");
		exit(-1);
    }
    double temp;
    for(int r=0;r<skip;r++){
        if(fscanf(fp, "%.17g", &temp) == 0){
            perror("Incorrect matrix dimensions");
            exit(-1);
        }
    }
    for (int i=0;i<M;i++){
        for(int j=0;j<N;j++){
            if(fscanf(fp, "%.17g", A + j + i*N) == 0){
                printf("Goes wrong on row %d and col %d\n", i, j);
			    perror("Incorrect matrix dimensions");
			    exit(-1);
		    }
        }
	}
    fclose(fp);
}

/**
 * @brief Function that gathers all parts of Q from other processes.
 * @param Q The array that will hold the result if the calling process is 0, otherwise the part of Q.
 * @param m The number of rows in the part.
 * @param n The number of columns in the part.
 * @param myid The rank of the calling process.
 * @param nprocs The number of processes.
 * @param start The start indices of all row blocks for all processes.
 * @param comm The MPI communicator.
 */
void GatherQ(double * Q, int m, int n, const int myid, const int nprocs, int * start, MPI_Comm comm){
  if(nprocs > 1){
    double * temp;
    if(myid>0){
        MPI_Send(&m, 1, MPI_INT, 0, 1, comm);
        MPI_Send(Q, m*n, MPI_DOUBLE, 0, 2, comm);
    }else{
      for(int p=1;p<nprocs;p++){
        MPI_Recv(&m, 1, MPI_INT, p, 1, comm, MPI_STATUS_IGNORE);
        temp = malloc(m*n*sizeof(double));
        MPI_Recv(temp, m*n, MPI_DOUBLE, p, 2, comm, MPI_STATUS_IGNORE);
        for(int j = 0;j<m;j++){
            for(int i=0;i<n;i++){
                Q[(start[p] + j)*n + i] = temp[i*m + j];
            }
        }
        free(temp);
      }
    }
  }
}
