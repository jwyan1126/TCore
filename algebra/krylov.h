#ifndef KRYLOV_H
#define KRYLOV_H
#include"vec.h"
#include"mat.h"

typedef int (*MATSOLVER)(VEC *x, const MAT *A, const VEC *b, int max_iter_num);

// Solve Ax=b with BiCGSTAB without preconditioner
// x		:initial guess values, and store the solution
// A		:matrix A
// b		:matrix b
// return	: 0 when converged
// 		: -1 when reach max iteration num
int bicgstab(VEC *x, const MAT *A, const VEC *b, int max_iter_num);

#endif
