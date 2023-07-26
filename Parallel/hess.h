/**
 * @file hess.h
 * @brief Header file related to the calculation of the upper Hessenberg matrix as part 
 * of CA-Arnoldi, part of Thesis project in High Performance Computing at Trinity College Dublin.
 * @author Hanne Meuwissen (meuwissh@tcd.ie)
 * @version 2.0
 * @date 2023-06-19
 */
#ifndef HESS_H_7BZWC1IU
#define HESS_H_7BZWC1IU

void calc_hess_on_transpose(double * H_, double * R_, double * B_, const int n, const int m);
void get_R(double * R, double * R_, const int n);
void set_B_(double *B_, const int s);
void update_hess_on_transpose(double ** H, double * mathcalR_, double * R_, const int s, const int k);
int breakdown_check(double *H, const int s, const int k, const double tol);

#endif /* end of include guard: HESS_H_7BZWC1IU */