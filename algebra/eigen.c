#include"eigen.h"
#include<stdlib.h>
#include"mvop.h"

void pow_iter(double *lambda, VEC *x, const MAT *A, int max_iter_num)
{
	#ifdef DEBUG
	if(x->size != A->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = x->size;
	VEC *x_last = vec_create(size);
	double val_last;
	double val = 1.0;
	vec_normalize(x);
	int counter = 0;
	while(counter < max_iter_num){
		counter++;
		vec_copy(x_last, x);
		val_last = val;
		mat_vec_mult_assign(x, A, x_last, 1.0/val_last);
		val = val_last * vec_inner_prod(x,x) / vec_inner_prod(x,x_last);
		double res = vec_res_2norm(x,x_last);
		//printf("Res=%g\n", res);
		if(res < 1e-5) break;
	}
	*lambda = val;
	vec_free(x_last);
	if(counter == max_iter_num){
		fprintf(stderr ,"pow_iter NOT converged.\n");
		exit(-1);
	}
}

void ipow_iter(double *lambda, VEC *x, const MAT *A, int max_iter_num)
{
	#ifdef DEBUG
	if(x->size != A->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = x->size;
	VEC *x_last = vec_create(size);
	double val_last;
	double val = 1.0;
	vec_normalize(x);
	int counter = 0;
	while(counter < max_iter_num){
		counter++;
		vec_copy(x_last, x);
		val_last = val;
		(*mat_solver)(x, A, x_last, 256);
		vec_scale(x, 1.0/val_last);
		val = val_last * vec_inner_prod(x,x) / vec_inner_prod(x,x_last);
		double res = vec_res_2norm(x,x_last);
		//printf("Res=%g\n", res);
		if(res < 1e-5) break;
	}
	*lambda = 1.0 / val;
	vec_free(x_last);
	if(counter == max_iter_num){
		fprintf(stderr, "ipow_iter NOT converged.\n");
		exit(-1);
	}
}

void spow_iter(double *lambda, VEC *x, const MAT *A, double sigma, int max_iter_num)
{
	#ifdef DEBUG
	if(x->size != A->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = x->size;
	VEC *x_last = vec_create(size);
	MAT *sA = mat_create(size);
	mat_copy(sA, A);
	for(size_t i=0; i<size; ++i)
		sA->vals[i][i] -= sigma;
	double val_last;
	double val = 1.0;
	vec_normalize(x);
	int counter = 0;
	while(counter < max_iter_num){
		counter++;
		vec_copy(x_last, x);
		val_last = val;
		(*mat_solver)(x, sA, x_last, 256);
		vec_scale(x, 1.0/val_last);
		val = val_last * vec_inner_prod(x,x) / vec_inner_prod(x,x_last);
		double res = vec_res_2norm(x,x_last);
		//printf("Res=%g\n", res);
		if(res < 1e-5) break;
	}
	*lambda = 1.0 / val + sigma;
	mat_free(sA);
	vec_free(x_last);
	if(counter == max_iter_num){
		fprintf(stderr, "spow_iter NOT converged.\n");
		exit(-1);
	}
}

void gipow_iter(double *lambda, VEC *x, const MAT *A, const MAT *B, int max_iter_num)
{
	#ifdef DEBUG
	if(x->size != A->size || x->size != B->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = x->size;
	VEC *x_last = vec_create(size);
	VEC *rvec = vec_create(size);
	double val_last;
	double val = 1.0;
	vec_normalize(x);
	int counter = 0;
	while(counter < max_iter_num){
		counter++;
		vec_copy(x_last, x);
		val_last = val;
		mat_vec_mult_assign(rvec, B, x_last, 1.0);
		(*mat_solver)(x, A, rvec, 256);
		vec_scale(x, 1.0/val_last);
		val = val_last * vec_inner_prod(x,x) / vec_inner_prod(x,x_last);
		double res = vec_res_2norm(x,x_last);
		//printf("Res=%g\n", res);
		if(res < 1e-5) break;
	}
	*lambda = 1.0 / val;
	vec_free(rvec);
	vec_free(x_last);
	if(counter == max_iter_num){
		fprintf(stderr, "gipow_iter NOT converged.\n");
		exit(-1);
	}
}

void gspow_iter(double *lambda, VEC *x, const MAT *A, const MAT *B, double sigma, int max_iter_num)
{
	#ifdef DEBUG
	if(x->size != A->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = x->size;
	VEC *x_last = vec_create(size);
	MAT *sA = mat_create(size);
	VEC *rvec = vec_create(size);
	mat_copy(sA, A);
	mat_adds(sA, B, -sigma);
	double val_last;
	double val = 1.0;
	vec_normalize(x);
	int counter = 0;
	while(counter < max_iter_num){
		counter++;
		vec_copy(x_last, x);
		val_last = val;
		mat_vec_mult_assign(rvec, B, x_last, 1.0);
		(*mat_solver)(x, sA, rvec, 256);
		vec_scale(x, 1.0/val_last);
		val = val_last * vec_inner_prod(x,x) / vec_inner_prod(x,x_last);
		double res = vec_res_2norm(x,x_last);
		//printf("Res=%g\n", res);
		if(res < 1e-5) break;
	}
	*lambda = 1.0 / val + sigma;
	vec_free(rvec);
	mat_free(sA);
	vec_free(x_last);
	if(counter == max_iter_num){
		fprintf(stderr, "gspow_iter NOT converged.\n");
		exit(-1);
	}
}
