#ifndef SSOL_H
#define SSOL_H

#include<stddef.h>

typedef struct
{
	double keff;
	size_t rt_size;
	size_t eg_size;
	double *flux;
} SSOL;

SSOL *ssol_create(size_t eg_size, size_t rt_size);

void ssol_free(SSOL *ssol);

#endif
