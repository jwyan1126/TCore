#include<stdio.h>
#include"vec.h"
#include"mat.h"
#include"mvop.h"
#include"ksp.h"
#include<time.h>
#include<stdlib.h>
#include"eigen.h"

MATSOLVER mat_solver;

int main()
{
	const size_t size = 1024;
	MAT *A = mat_create(size);
	mat_set_rand(A);
	VEC *b = vec_create(size);
	vec_set_ones(b);
	VEC *x = vec_create(size);
	LU_solve(x, A, b);
	VEC *bp = vec_create(size);
	mat_vec_mult_assign(bp, A, x, 1.0);
	double res = vec_res_2norm(b,bp);
	printf("res=%g\n", res);
	vec_free(bp);
	vec_free(x);
	vec_free(b);
	mat_free(A);
	return 0;
}
