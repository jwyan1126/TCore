#ifndef KSP_H
#define KSP_H
#include"vec.h"
#include"mat.h"

// x		:initial guess values, and store the solution
// A		:matrix A
// b		:matrix b
// return	: 0 when converged
// 		: -1 when reach max iteration num
typedef void (*MATSOLVER)(VEC *x, const MAT *A, const VEC *b, int max_iter_num);

// Solve Ax=b with BiCGSTAB without preconditioner
void bicgstab(VEC *x, const MAT *A, const VEC *b, int max_iter_num);

// Solve Ax=b with gauss_seidel
void gauss_seidel(VEC *x, const MAT *A, const VEC *b, int max_iter_num);

void LU_solve(VEC *x, const MAT *A, const VEC *b);

void LU_decomposition(MAT *LU, const MAT *A);

#endif
