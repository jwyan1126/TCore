#ifndef KSP_H
#define KSP_H
#include"vec.h"
#include"mat.h"

// x		:initial guess values, and store the solution
// A		:matrix A
// b		:matrix b
// return	: 0 when converged
// 		: -1 when reach max iteration num
typedef int (*MATSOLVER)(VEC *x, const MAT *A, const VEC *b, int max_iter_num);

// Solve Ax=b with BiCGSTAB without preconditioner
int bicgstab(VEC *x, const MAT *A, const VEC *b, int max_iter_num);

// Solve Ax=b with gauss_seidel
int gauss_seidel(VEC *x, const MAT *A, const VEC *b, int max_iter_num);

#endif
