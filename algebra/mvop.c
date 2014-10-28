#include"mvop.h"
#include"vec.h"
#include"mat.h"
#include<stdlib.h>

void mat_vec_mult_assign(VEC *tar_vec, const MAT *mat, const VEC *vec, double alpha)
{
	#ifdef DEBUG
	if(tar_vec->size != vec->size || tar_vec->size != mat->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = tar_vec->size;
	for(size_t i=0; i<size; ++i){
		tar_vec->vals[i] = 0.0;
		for(size_t j=0; j<size; ++j)
			tar_vec->vals[i] += alpha * mat->vals[i][j] * vec->vals[j];
	}
}

void mat_vec_Tmult_assign(VEC *tar_vec, const MAT *mat, const VEC *vec, double alpha)
{
	#ifdef DEBUG
	if(tar_vec->size != vec->size || tar_vec->size != mat->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = tar_vec->size;
	for(size_t i=0; i<size; ++i){
		tar_vec->vals[i] = 0.0;
		for(size_t j=0; j<size; ++j)
			tar_vec->vals[i] += alpha * mat->vals[j][i] * vec->vals[j];
	}
}

void mat_vec_mult_update(VEC *tar_vec, const MAT *mat, const VEC *vec, double alpha)
{
	#ifdef DEBUG
	if(tar_vec->size != vec->size || tar_vec->size != mat->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = tar_vec->size;
	for(size_t i=0; i<size; ++i){
		for(size_t j=0; j<size; ++j)
			tar_vec->vals[i] += alpha * mat->vals[i][j] * vec->vals[j];
	}
}

void mat_vec_Tmult_update(VEC *tar_vec, const MAT *mat, const VEC *vec, double alpha)
{
	#ifdef DEBUG
	if(tar_vec->size != vec->size || tar_vec->size != mat->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = tar_vec->size;
	for(size_t i=0; i<size; ++i){
		for(size_t j=0; j<size; ++j)
			tar_vec->vals[i] += alpha * mat->vals[j][i] * vec->vals[j];
	}
}
