/**
 * @file main.c
 * @brief Main function for thesis.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 1.0
 * @date 2023-05-26
 */

#include<stdlib.h>                                                             
#include<stdio.h>
#include<string.h> 
#include"graph.h"
#include"sparse.h"                                                       
                                                                                
int main(void){
    /* Check graph sparse transition matrix generation & spmv */
    // int n = 5;
    // int nnz = 3;                                                                        
    // sparse_CSR test = generate_regular_graph_trans_csr(n, nnz);
    // print_CSR(&test);
    // double * v = malloc(n*sizeof(double));
    // for(int i=0;i<n;i++){v[i] = 0.5;}
    // double * result = malloc(n*sizeof(double));
    // spmv(test, v, n, result);
    // print_vector(result, n);
    // free(test.rowptrs);
    // free(test.colindex);
    // free(test.values);

    /* Check serial arnoldi */
    int n = 3;
    int nnz_per_row = 2;
    sparse_CSR T;
    T.nrows = n;
    T.ncols = n;
    T.nnz = n*nnz_per_row;
    T.rowptrs = malloc((n+1)*sizeof(int));
    T.colindex = malloc(T.nnz*sizeof(int));
    T.values = malloc(T.nnz*sizeof(double));
    for(int i=0;i<4;i++){
        T.rowptrs[i] = i*2;
    }
    T.colindex[0] = 0;
    T.colindex[1] = 1;
    T.colindex[2] = 1;
    T.colindex[3] = 2;
    T.colindex[4] = 0;
    T.colindex[5] = 1;
    T.values[0] = 1.0;
    T.values[1] = 2.0;
    T.values[2] = 1.0;
    T.values[3] = 2.0;
    T.values[4] = 2.0;
    T.values[5] = 1.0;
    // print_CSR(&T);
    int m = 3;
    double * Q = malloc(n*m*sizeof(double));
    double *b = malloc(n*sizeof(double));
    for(int i=0;i<3;i++){
        b[i] = 1.0/sqrt(3.0);
    }
    Arnoldi(T, b, n, Q, m);
    print_matrix(Q, n, m);

}     