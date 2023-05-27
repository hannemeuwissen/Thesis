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
    sparse_CSR test = generate_regular_graph_trans_csr(4, 1);
    print_CSR(&test);
    double * v = malloc(4*sizeof(double));
    for(int i=0;i<4;i++){v[i] = 1.0;}
    double * result = malloc(4*sizeof(double));
    memset(result, 0, 4*sizeof(double));
    spmv(test, v, result);
    print_vector(result);
    free(test.rowptrs);
    free(test.colindex);
    free(test.values);
}     