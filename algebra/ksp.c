#include"ksp.h"
#include<stdlib.h>
#include"mvop.h"
#include<math.h>
#include<stddef.h>

void bicgstab(VEC *x, const MAT *A, const VEC *b, int max_iter_num)
{
	#ifdef DEBUG
	if(x->size != A->size || x->size != b->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = x->size;
	VEC *r = vec_create(size);
	VEC *r_hat = vec_create(size);
	VEC *v = vec_create(size);
	VEC *p = vec_create(size);
	VEC *s = vec_create(size);
	VEC *t = vec_create(size);
	VEC *x_last = vec_create(size);
	vec_copy(r, b);
	mat_vec_mult_update(r, A, x, -1.0);
	vec_copy(r_hat, r);
	double rou = 1.0;
	double rou_last;
	double alpha = 1.0;
	double omega = 1.0;
	int counter = 0;
	while(counter < max_iter_num){
		counter++;
		rou_last = rou;
		rou = vec_inner_prod(r_hat, r);
		double beta = (rou/rou_last)*(alpha/omega);
		vec_adds_assign(s, p, v, 1.0, -omega);
		vec_adds_assign(p, r, s, 1.0, beta);
		mat_vec_mult_assign(v, A, p, 1.0);
		alpha = rou / vec_inner_prod(r_hat, v);
		vec_adds_assign(s, r, v, 1.0, -alpha);
		mat_vec_mult_assign(t, A, s, 1.0);
		omega = vec_inner_prod(t,s) / vec_inner_prod(t,t);
		vec_copy(x_last, x);
		vec_adds(x, p, alpha);
		vec_adds(x, s, omega);
		double res = vec_res_2norm(x, x_last);
		if(res < 1e-6) break;
		vec_adds_assign(r, s, t, 1.0, -omega);
	}
	vec_free(x_last);
	vec_free(t);
	vec_free(s);
	vec_free(p);
	vec_free(v);
	vec_free(r_hat);
	vec_free(r);
	if(counter == max_iter_num){
		fprintf(stderr, "BICGSTAB NOT converged.\n");
		exit(-1);
	}
}

void gauss_seidel(VEC *x, const MAT *A, const VEC *b, int max_iter_num)
{
	#ifdef DEBUG
	if(x->size != A->size || x->size != b->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = x->size;
	int counter = 0;
	VEC *x_last = vec_create(size);
	while(counter < max_iter_num){
		counter++;
		vec_copy(x_last, x);
		for(size_t i=0; i<size; ++i){
			double sum = 0.0;
			for(size_t j=0; j<i; ++j)
				sum += mat_get(A, i, j)*vec_get(x,j);
			for(size_t j=i+1; j<size; ++j)
				sum += mat_get(A, i, j)*vec_get(x,j);
			vec_set(x, i, (vec_get(b,i) - sum) / mat_get(A,i,i));
		}
		double res = vec_res_2norm(x, x_last);
		//printf("%g\n", res);
		if(res < 1e-6) break;
	}
	vec_free(x_last);
	if(counter == max_iter_num){
		fprintf(stderr, "GS NOT converged.\n");
		exit(-1);
	}
}

void LU_solve(VEC *x, const MAT *A, const VEC *b)
{
	#ifdef DEBUG
	if(x->size != A->size || x->size != b->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = x->size;
	MAT *LU = mat_create(size);
	LU_decomposition(LU, A);
	for(size_t i=0; i<size; ++i){
		vec_set(x, i, vec_get(b, i));
		for(size_t j=0; j<i; ++j)
			vec_set(x, i, vec_get(x,i) - mat_get(LU, i, j)*vec_get(x, j));
	}
	for(size_t i=size; i>0; --i){
		for(size_t j=i; j<size; ++j)
			vec_set(x, i-1, vec_get(x,i-1) - mat_get(LU, i-1, j)*vec_get(x, j));
		vec_set(x, i-1, vec_get(x, i-1) / mat_get(LU, i-1, i-1));
	}
	free(LU);
}

void LU_decomposition(MAT *LU, const MAT *A)
{
	#ifdef DEBUG
	if(LU->size != A->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = A->size;
	// init
	for(size_t j=0; j<size; ++j)
		mat_set(LU, 0, j, mat_get(A, 0, j));
	for(size_t i=1; i<size; ++i){
		double pivot = mat_get(LU, 0, 0);
		#ifdef DEBUG
		if(fabs(pivot) < 1e-6){
			fprintf(stderr, "LU failed.\n");
			exit(-1);
		}
		#endif
		mat_set(LU, i, 0, mat_get(A, i, 0)/pivot);
	}

	for(size_t p=1; p<size; ++p){
		// u
		for(size_t j=p; j<size; ++j){
			mat_set(LU, p, j, mat_get(A, p, j));
			for(size_t r=0; r<p; ++r)
				mat_set(LU, p, j, mat_get(LU, p, j) - mat_get(LU, p, r)*mat_get(LU, r, j));
		}
		// l
		for(size_t i=p+1; i<size; ++i){
			double pivot = mat_get(LU, p, p);
			#ifdef DEBUG
			if(fabs(pivot) < 1e-6){
				fprintf(stderr, "LU failed.\n");
				exit(-1);
			}
			#endif
			mat_set(LU, i, p, mat_get(A, i, p));
			for(size_t r=0; r<p; ++r)
				mat_set(LU, i, p, mat_get(LU, i, p) - mat_get(LU, i, r)*mat_get(LU, r, p));
			mat_set(LU, i, p, mat_get(LU, i, p) / pivot);
		}
	}
}
