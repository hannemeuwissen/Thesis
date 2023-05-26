/**
 * @file main.c
 * @brief Main function for thesis.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 1.0
 * @date 2023-05-26
 */

#include<stdlib.h>                                                             
#include<stdio.h>    
#include"graph.h"                                                          
                                                                                
int main(void){                                                                               
    sparse_CSR test = generate_regular_graph_trans_csr(5, 2);
    print_CSR(&test);
    free(test.rowptrs);
    free(test.colindex);
    free(test.values);
}     