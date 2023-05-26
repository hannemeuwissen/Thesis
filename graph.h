/**
 * @file graph.h
 * @brief Main function for thesis.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 1.0
 * @date 2023-05-26
 */
#ifndef GRAPH_H_TWYFUKH2
#define GRAPH_H_TWYFUKH2

typedef struct sparse_CSR {
    size_t nrows;
    size_t ncols;
    size_t nnz;
    size_t * rowptrs;
    size_t * colindex;
    double * values;
} sparse_CSR;

sparse_CSR generate_regular_graph_trans_csr(const size_t n, const size_t nnz_per_row);
void print_CSR(sparse_CSR * M);

#endif /* end of include guard: GRAPH_H_TWYFUKH2 */
