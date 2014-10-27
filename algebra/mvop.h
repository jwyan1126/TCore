#ifndef MVOP_H
#define MVOP_H

#include"vec.h"
#include"mat.h"

// tar_vec = tar_vec + alpha * mat * vec
void mat_vec_mult_update(VEC *tar_vec, const MAT *mat, const VEC *vec, double alpha);

// tar_vec = alpha * mat * vec
void mat_vec_mult_assign(VEC *tar_vec, const MAT *mat, const VEC *vec, double alpha);

// tar_vec = tar_vec + alpha * mat' * vec
void mat_vec_Tmult_update(VEC *tar_vec, const MAT *mat, const VEC *vec, double alpha);

// tar_vec = alpha * mat' * vec
void mat_vec_Tmult_assign(VEC *tar_vec, const MAT *mat, const VEC *vec, double alpha);

#endif
