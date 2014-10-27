#ifndef VEC_H
#define VEC_H

#include<stddef.h>
#include<stdio.h>

typedef struct
{
	size_t size;
	double *vals;
} VEC;

VEC *vec_create(size_t size);

VEC *vec_ref_create(size_t size, double *vals);

void vec_free(VEC *vec);

void vec_ref_free(VEC *vec);

double vec_get(const VEC *vec, size_t i);

void vec_set(VEC *vec, size_t i, double val);

void vec_set_zeros(VEC *vec);

void vec_set_ones(VEC *vec);

void vec_set_rand(VEC *vec);

void vec_copy(VEC *tar_vec, const VEC *src_vec);

void vec_printf(const VEC *vec);

void vec_fprintf(const VEC *vec, FILE *stream);

double vec_1norm(const VEC *vec);

double vec_2norm(const VEC *vec);

double vec_infnorm(const VEC *vec);

// 1norm(vec1 - vec2)
double vec_res_1norm(const VEC *vec1, const VEC *vec2);

// 2norm(vec1 - vec2)
double vec_res_2norm(const VEC *vec1, const VEC *vec2);

// infnorm(vec1 - vec2)
double vec_res_infnorm(const VEC *vec1, const VEC *vec2);

void vec_normalize(VEC *vec);

double vec_inner_prod(const VEC *vec1, const VEC *vec2);

void vec_abs(VEC *vec);

void vec_abs_assign(VEC *tar_vec, const VEC *vec);

void vec_scale(VEC *vec, double scale);

void vec_scale_assign(VEC *tar_vec, const VEC *vec, double scale);

void vec_add(VEC *vec, const VEC *from_vec);

void vec_add_assign(VEC *tar_vec, const VEC *vec1, const VEC *vec2);

// vec = vec + from_vec * scale
void vec_adds(VEC *vec, const VEC *from_vec, double scale);

// tar_vec = vec1*scale1 + vec2*scale2
void vec_adds_assign(VEC *tar_vec, const VEC *vec1, const VEC *vec2, double scale1, double scale2);

#endif
