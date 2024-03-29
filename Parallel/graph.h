/**
 * @file graph.h
 * @brief Header file related to the generation of regular graph transition matrices in sparse 
 * format to investigate CA-Arnoldi, part of Thesis project in High Performance Computing 
 * at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 3.0
 * @date 2023-05-26
 */
#ifndef GRAPH_H_TWYFUKH2
#define GRAPH_H_TWYFUKH2

#include"sparse.h"
sparse_CSR generate_regular_graph_part_csr(const int n, const int M, const int nnz_per_row, const int random);
sparse_CSR generate_irregular_graph_part_csr(const int n, const int M, const int min_nnz_per_row, const int max_nnz_per_row, const int random_ind);
sparse_CSR generate_irregular_csr(const int M, const int min_nnz, const int max_nnz);

#endif /* end of include guard: GRAPH_H_TWYFUKH2 */
