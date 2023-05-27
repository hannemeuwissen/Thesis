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
    int nrows;
    int ncols;
    int nnz;
    int * rowptrs;
    int * colindex;
    double * values;
} sparse_CSR;

sparse_CSR generate_regular_graph_trans_csr(const int n, const int nnz_per_row);
void print_CSR(sparse_CSR * M);
void spmv(sparse_CSR M, double * v, double * result);
void print_vector(double * v, const int len);

#endif /* end of include guard: GRAPH_H_TWYFUKH2 */
