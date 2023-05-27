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
                                                                                
int main(void){
    int n = 5;
    int nnz = 3;                                                                        
    sparse_CSR test = generate_regular_graph_trans_csr(n, nnz);
    print_CSR(&test);
    double * v = malloc(n*sizeof(double));
    for(int i=0;i<n;i++){v[i] = 0.5;}
    double * result = malloc(n*sizeof(double));
    spmv(test, v, n, result);
    print_vector(result, n);
    free(test.rowptrs);
    free(test.colindex);
    free(test.values);
}     