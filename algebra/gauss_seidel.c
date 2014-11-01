#include"ksp.h"
#include<stdlib.h>

int gauss_seidel(VEC *x, const MAT *A, const VEC *b, int max_iter_num)
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
	return counter == max_iter_num? -1 : 0;
}
