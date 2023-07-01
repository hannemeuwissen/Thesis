/**
 * @file parse.h
 * @brief Header for thesis at Trinity College Dublin.
 * @author Hanne Meuwissen (22307813)
 * @version 1.0
 * @date 2023-06-02
 */
#ifndef PARSE_H_7BZWC1IU
#define PARSE_H_7BZWC1IU

void parse_command_line_regular(const int argc, char * const *argv, int * M, int * nnz, char * filename_v, int * degree, int * s, MPI_Comm comm);
void parse_command_line_irregular(const int argc, char * const *argv, char * filename_A, int * M, char * filename_v, int * degree, int * s, MPI_Comm comm);

#endif /* end of include guard: PARSE_H_7BZWC1IU */