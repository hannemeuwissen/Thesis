/**
 * @file hess.h
 * @brief Header for thesis at Trinity College Dublin.
 * @author Hanne Meuwissen (22307813)
 * @version 1.0
 * @date 2023-06-19
 */
#ifndef HESS_H_7BZWC1IU
#define HESS_H_7BZWC1IU

void calc_hess(double * H_, double * R_, double * B_, double * R, const int n, const int m);
void get_R(double * R, double * R_, const int n);

#endif /* end of include guard: HESS_H_7BZWC1IU */