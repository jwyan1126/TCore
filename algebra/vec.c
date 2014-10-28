#include"vec.h"
#include<stdlib.h>
#include<math.h>

VEC *vec_create(size_t size)
{
	VEC *vec = malloc(sizeof(VEC));
	vec->size = size;
	vec->vals = calloc(size, sizeof(double));
	return vec;
}

VEC *vec_ref_create(size_t size, double *vals)
{
	VEC *vec = malloc(sizeof(VEC));
	vec->size = size;
	vec->vals = vals;
	return vec;
}

void vec_free(VEC *vec)
{
	free(vec->vals);
	free(vec);
}

void vec_ref_free(VEC *vec)
{
	free(vec);
}

inline double vec_get(const VEC *vec, size_t i)
{
	#ifdef DEBUG
	if(i >= vec->size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return vec->vals[i];
}

inline void vec_set(VEC *vec, size_t i, double val)
{
	#ifdef DEBUG
	if(i >= vec->size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	vec->vals[i] = val;
}

void vec_set_zeros(VEC *vec)
{
	for(size_t i=0; i<vec->size; ++i)
		vec->vals[i] = 0.0;
}

void vec_set_ones(VEC *vec)
{
	for(size_t i=0; i<vec->size; ++i)
		vec->vals[i] = 1.0;
}

void vec_copy(VEC *tar_vec, const VEC *src_vec)
{
	#ifdef DEBUG
	if(tar_vec->size != src_vec->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	for(size_t i=0; i<tar_vec->size; ++i)
		tar_vec->vals[i] = src_vec->vals[i];
}

void vec_printf(const VEC *vec)
{
	vec_fprintf(vec,stdout);
}

void vec_fprintf(const VEC *vec, FILE *stream)
{
	for(size_t i=0; i<vec->size; ++i)
		fprintf(stream, "%g\n", vec_get(vec,i));
	fprintf(stream, "\n");
}

void vec_set_rand(VEC *vec)
{
	for(size_t i=0; i<vec->size; ++i)
		vec->vals[i] = rand() / (RAND_MAX + 1.0);
}

double vec_1norm(const VEC *vec)
{
	size_t size = vec->size;
	double sum = 0.0;
	for(size_t i=0; i<size; ++i)
		sum += fabs(vec->vals[i]);
	return sum;
}

double vec_2norm(const VEC *vec)
{
	size_t size = vec->size;
	double sum = 0.0;
	for(size_t i=0; i<size; ++i)
		sum += vec->vals[i]*vec->vals[i];
	return sqrt(sum);
}

double vec_infnorm(const VEC *vec)
{
	double max = 0.0;
	size_t size = vec->size;
	for(size_t i=0; i<size; ++i){
		double max = 0.0;
		double val = fabs(vec->vals[i]);
		if(val > max){
			max = val;
		}
	}
	return max;
}

double vec_res_1norm(const VEC *vec1, const VEC *vec2)
{
	#ifdef DEBUG
	if(vec1->size != vec2->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = vec1->size;
	double sum = 0.0;
	for(size_t i=0; i<size; ++i)
		sum += fabs(vec1->vals[i] - vec2->vals[i]);
	return sum;
}

double vec_res_2norm(const VEC *vec1, const VEC *vec2)
{
	#ifdef DEBUG
	if(vec1->size != vec2->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = vec1->size;
	double sum = 0.0;
	for(size_t i=0; i<size; ++i){
		double val = vec1->vals[i] - vec2->vals[i];
		sum += val*val;
	}
	return sqrt(sum);
}

double vec_res_infnorm(const VEC *vec1, const VEC *vec2)
{
	#ifdef DEBUG
	if(vec1->size != vec2->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = vec1->size;
	double max = 0.0;
	for(size_t i=0; i<size; ++i){
		double val = fabs(vec1->vals[i] - vec2->vals[i]);
		if(val > max){
			max = val;
		}
		
	}
	return max;
}

void vec_normalize(VEC *vec)
{
	double len = vec_2norm(vec);
	vec_scale(vec, 1.0/len);
}

double vec_inner_prod(const VEC *vec1, const VEC *vec2)
{
	#ifdef DEBUG
	if(vec1->size != vec2->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = vec1->size;
	double sum = 0.0;
	for(size_t i=0; i<size; ++i)
		sum += vec1->vals[i] * vec2->vals[i];
	return sum;
}

void vec_abs(VEC *vec)
{
	size_t size = vec->size;
	for(size_t i=0; i<size; ++i)
		vec->vals[i] = fabs(vec->vals[i]);
}

void vec_abs_assign(VEC *tar_vec, const VEC *vec)
{
	#ifdef DEBUG
	if(tar_vec->size != vec->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = tar_vec->size;
	for(size_t i=0; i<size; ++i)
		tar_vec->vals[i] = fabs(vec->vals[i]);
}

void vec_scale(VEC *vec, double scale)
{
	size_t size = vec->size;
	for(size_t i=0; i<size; ++i)
		vec->vals[i] *= scale;
}

void vec_scale_assign(VEC *tar_vec, const VEC *vec, double scale)
{
	#ifdef DEBUG
	if(tar_vec->size != vec->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = tar_vec->size;
	for(size_t i=0; i<size; ++i)
		tar_vec->vals[i] = vec->vals[i] * scale;
}

void vec_add(VEC *vec, const VEC *from_vec)
{
	#ifdef DEBUG
	if(vec->size != from_vec->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = vec->size;
	for(size_t i=0; i<size; ++i)
		vec->vals[i] += from_vec->vals[i];
}

void vec_add_assign(VEC *tar_vec, const VEC *vec1, const VEC *vec2)
{
	#ifdef DEBUG
	if(tar_vec->size != vec1->size || tar_vec->size != vec2->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = tar_vec->size;
	for(size_t i=0; i<size; ++i)
		tar_vec->vals[i] = vec1->vals[i] + vec2->vals[i];
}

void vec_adds(VEC *vec, const VEC *from_vec, double scale)
{
	#ifdef DEBUG
	if(vec->size != from_vec->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = vec->size;
	for(size_t i=0; i<size; ++i)
		vec->vals[i] += scale * from_vec->vals[i];
}

void vec_adds_assign(VEC *tar_vec, const VEC *vec1, const VEC *vec2, double scale1, double scale2)
{
	#ifdef DEBUG
	if(tar_vec->size != vec1->size || tar_vec->size != vec2->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = tar_vec->size;
	for(size_t i=0; i<size; ++i)
		tar_vec->vals[i] = scale1*vec1->vals[i] + scale2*vec2->vals[i];
}
