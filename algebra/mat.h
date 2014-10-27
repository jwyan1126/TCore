#ifndef MAT_H
#define MAT_H

#include<stddef.h>
#include<stdio.h>

// Square matrix with dense storage
typedef struct
{
	size_t size;
	double **vals;
} MAT;

MAT *mat_create(size_t size);

MAT *mat_ref_create(size_t size, double **vals);

void mat_free(MAT *mat);

void mat_ref_free(MAT *mat);

double mat_get(const MAT *mat, size_t i, size_t j);

void mat_set(MAT *mat, size_t i, size_t j, double val);

void mat_set_zeros(MAT *mat);

void mat_set_ones(MAT *mat);

void mat_set_identity(MAT *mat);

void mat_set_rand(MAT *mat);

void mat_copy(MAT *tar_mat, const MAT *src_mat);

void mat_transpose_assign(MAT *tar_mat, const MAT *mat);

void mat_printf(const MAT *mat);

void mat_fprintf(const MAT *mat, FILE *stream);

void mat_abs(MAT *mat);

void mat_abs_assign(MAT *tar_mat, const MAT *mat);

void mat_scale(MAT *mat, double scale);

void mat_scale_assign(MAT *tar_mat, const MAT *mat, double scale);

void mat_add(MAT *mat, const MAT *from_mat);

void mat_add_assign(MAT *tar_mat, const MAT *mat1, const MAT *mat2);

// mat = mat + from_mat * scale
void mat_adds(MAT *mat, const MAT *from_mat, double scale);

// tar_mat = mat1*scale1 + mat2*scale2
void mat_adds_assign(MAT *tar_mat, const MAT *mat1, const MAT *mat2, double scale1, double scale2);

#endif
