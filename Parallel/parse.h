/**
 * @file parse.h
 * @brief Header file related to parsing the input parameters to investigate CA-Arnoldi, 
 * part of Thesis project in High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 2.0
 * @date 2023-06-02
 */
#ifndef PARSE_H_7BZWC1IU
#define PARSE_H_7BZWC1IU

void parse_command_line_regular(const int argc, char * const *argv, int * M, int * nnz, char * filename_v, int * degree, int * s, MPI_Comm comm);
void parse_command_line_irregular(const int argc, char * const *argv, int * M, int * min_nnz, int * max_nnz, char * filename_v, int * degree, int * s, MPI_Comm comm);
void parse_command_line_lb(const int argc, char * const *argv, char * filename_A, char * filename_v, int * degree, int * s, int * lb, int * q, int * h, MPI_Comm comm);
#endif /* end of include guard: PARSE_H_7BZWC1IU */