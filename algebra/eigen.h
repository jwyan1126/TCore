#ifndef EIGEN_H
#define EIGEN_H
#include"mat.h"
#include"vec.h"
#include<stddef.h>
#include"ksp.h"

extern MATSOLVER mat_solver;

// A*x = lambda*x
// Conventional power iteration
// lambda	:output eigenvalue
// x		:initial guess values, and store the solution
// A		:matrix A
// return	: 0 when converged
// 		: -1 when reach max iteration num
int pow_iter(double *lambda, VEC *x, const MAT *A, int max_iter_num);

// A*x = lambda*x
// Inverse power iteration
int ipow_iter(double *lambda, VEC *x, const MAT *A, int max_iter_num);

// A*x = lambda*x
// Shifted power iteration
int spow_iter(double *lambda, VEC *x, const MAT *A, double sigma, int max_iter_num);

// A*x = lambda*B*x
// Generalized inverse power iteration
int gipow_iter(double *lambda, VEC *x, const MAT *A, const MAT *B, int max_iter_num);

// A*x = lambda*B*x
// Generalized shifted power iteration
int gspow_iter(double *lambda, VEC *x, const MAT *A, const MAT *B, double sigma, int max_iter_num);

#endif
