#include"krylov.h"
#include<stdlib.h>
#include"mvop.h"
#include<stddef.h>

int bicgstab(VEC *x, const MAT *A, const VEC *b, int max_iter_num)
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
	return counter == max_iter_num? -1 : 0;
}
