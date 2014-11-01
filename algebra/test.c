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
	mat_solver = bicgstab;
	//srand((int)time(NULL));
	size_t size = 16;
	MAT *A = mat_create(size);
	MAT *B = mat_create(size);
	mat_set_rand(A);
	mat_set_rand(B);
	MAT *tA = mat_create(size);
	mat_transpose_assign(tA,A);
	mat_add(A,tA);
	mat_transpose_assign(tA,B);
	mat_add(B,tA);
	mat_free(tA);
	for(size_t i=0; i<size; ++i){
		mat_set(A,i,i,size);
		mat_set(B,i,i,size/2);
	}
	double lambda;
	VEC *x = vec_create(size);
	vec_set_ones(x);
	if(gipow_iter(&lambda, x, A, B, 256))
		printf("Not converged.\n");
	else
		printf("Converged.\n");
	printf("eigenvalue=%g\n", lambda);
	vec_free(x);
	mat_free(B);
	mat_free(A);
	return 0;
}
