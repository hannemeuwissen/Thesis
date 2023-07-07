/**
 * @file graph.h
 * @brief Header file related to generating sparse transition matrices of random graphs, 
 * part of Thesis project in High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 1.0
 * @date 2023-05-26
 */
#ifndef GRAPH_H_TWYFUKH2
#define GRAPH_H_TWYFUKH2

#include"sparse.h"
sparse_CSR generate_regular_graph_trans_csr(const int n, const int nnz_per_row);

#endif /* end of include guard: GRAPH_H_TWYFUKH2 */
