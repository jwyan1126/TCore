#ifndef CAL_BGK_H
#define CAL_BGK_H

#include<stddef.h>

int delta_func(size_t g, size_t from_g);

double cal_bgk(size_t g, size_t from_g, double *Dgk, double *srgk, double *chigk, double **ssgk, double *vsfgk, double keff);

#endif
