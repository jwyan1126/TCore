#include"mat.h"

#include<math.h>
#include<stdlib.h>

MAT *mat_create(size_t size)
{
	MAT *mat = malloc(sizeof(MAT));
	mat->size = size;
	mat->vals = malloc(size * sizeof(double *));
	for(size_t i=0; i<size; ++i)
		mat->vals[i] = calloc(size,sizeof(double));
	return mat;
}

MAT *mat_ref_create(size_t size, double **vals)
{
	MAT *mat = malloc(sizeof(MAT));
	mat->size = size;
	mat->vals = vals;
	return mat;
}

void mat_free(MAT *mat)
{
	for(size_t i=0; i<mat->size; ++i)
		free(mat->vals[i]);
	free(mat->vals);
	free(mat);
}

void mat_ref_free(MAT *mat)
{
	free(mat);
}

double mat_get(const MAT *mat, size_t i, size_t j)
{
	#ifdef DEBUG
	if(i >= mat->size || j >= mat->size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return mat->vals[i][j];
}

void mat_set(MAT *mat, size_t i, size_t j, double val)
{
	#ifdef DEBUG
	if(i >= mat->size || j >= mat->size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	mat->vals[i][j] = val;
}

void mat_set_zeros(MAT *mat)
{
	size_t size = mat->size;
	for(size_t i=0; i<size; ++i)
		for(size_t j=0; j<size; ++j)
			mat->vals[i][j] = 0.0;
}

void mat_set_ones(MAT *mat)
{
	size_t size = mat->size;
	for(size_t i=0; i<size; ++i)
		for(size_t j=0; j<size; ++j)
			mat->vals[i][j] = 1.0;
}

void mat_set_identity(MAT *mat)
{
	size_t size = mat->size;
	for(size_t i=0; i<size; ++i)
		for(size_t j=0; j<size; ++j)
			mat->vals[i][j] = (i == j) ? 1.0 : 0.0;
}

void mat_set_rand(MAT *mat)
{
	size_t size = mat->size;
	for(size_t i=0; i<size; ++i)
		for(size_t j=0; j<size; ++j)
			mat->vals[i][j] = rand() / (RAND_MAX + 1.0);
}

void mat_copy(MAT *tar_mat, const MAT *src_mat)
{
	#ifdef DEBUG
	if(tar_mat->size != src_mat->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = tar_mat->size;
	for(size_t i=0; i<size; ++i)
		for(size_t j=0; j<size; ++j)
			tar_mat->vals[i][j] = src_mat->vals[i][j];
}

void mat_transpose_assign(MAT *tar_mat, const MAT *mat)
{
	#ifdef DEBUG
	if(tar_mat->size != mat->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = tar_mat->size;
	for(size_t i=0; i<size; ++i)
		for(size_t j=0; j<size; ++j)
			tar_mat->vals[j][i] = mat->vals[i][j];
}

void mat_printf(const MAT *mat)
{
	mat_fprintf(mat, stdout);
}

void mat_fprintf(const MAT *mat, FILE *stream)
{
	for(size_t i=0; i<mat->size; ++i){
		for(size_t j=0; j<mat->size; ++j){
			fprintf(stream, "%g\t", mat->vals[i][j]);
		}
		fprintf(stream, "\n");
	}
	fprintf(stream, "\n");
}

void mat_abs(MAT *mat)
{
	size_t size = mat->size;
	for(size_t i=0; i<size; ++i)
		for(size_t j=0; j<size; ++j)
			mat->vals[i][j] = fabs(mat->vals[i][j]);
}

void mat_abs_assign(MAT *tar_mat, const MAT *mat)
{
	#ifdef DEBUG
	if(tar_mat->size != mat->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = mat->size;
	for(size_t i=0; i<size; ++i)
		for(size_t j=0; j<size; ++j)
			tar_mat->vals[i][j] = fabs(mat->vals[i][j]);
}

void mat_scale(MAT *mat, double scale)
{
	size_t size = mat->size;
	for(size_t i=0; i<size; ++i)
		for(size_t j=0; j<size; ++j)
			mat->vals[i][j] *= scale;
}

void mat_scale_assign(MAT *tar_mat, const MAT *mat, double scale)
{
	#ifdef DEBUG
	if(tar_mat->size != mat->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = mat->size;
	for(size_t i=0; i<size; ++i)
		for(size_t j=0; j<size; ++j)
				tar_mat->vals[i][j] = mat->vals[i][j] * scale;
}

void mat_add(MAT *mat, const MAT *from_mat)
{
	#ifdef DEBUG
	if(mat->size != from_mat->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = mat->size;
	for(size_t i=0; i<size; ++i)
		for(size_t j=0; j<size; ++j)
			mat->vals[i][j] += from_mat->vals[i][j];
}

void mat_add_assign(MAT *tar_mat, const MAT *mat1, const MAT *mat2)
{
	#ifdef DEBUG
	if(tar_mat->size != mat1->size || tar_mat->size != mat2->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = tar_mat->size;
	for(size_t i=0; i<size; ++i)
		for(size_t j=0; j<size; ++j)
			tar_mat->vals[i][j] = mat1->vals[i][j] + mat2->vals[i][j];
}

void mat_adds(MAT *mat, const MAT *from_mat, double scale)
{
	#ifdef DEBUG
	if(mat->size != from_mat->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = mat->size;
	for(size_t i=0; i<size; ++i)
		for(size_t j=0; j<size; ++j)
			mat->vals[i][j] += scale * from_mat->vals[i][j];
}

void mat_adds_assign(MAT *tar_mat, const MAT *mat1, const MAT *mat2, double scale1, double scale2)
{
	#ifdef DEBUG
	if(tar_mat->size != mat1->size || tar_mat->size != mat2->size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t size = tar_mat->size;
	for(size_t i=0; i<size; ++i)
		for(size_t j=0; j<size; ++j)
			tar_mat->vals[i][j] = scale1 * mat1->vals[i][j] + scale2 * mat2->vals[i][j];
}
