#ifndef SSOL_H
#define SSOL_H

#include<stddef.h>
#include"flux.h"

typedef struct
{
	size_t eg_size;
	size_t rt_size;
	size_t xm_size;
	size_t ym_size;
	size_t zm_size;
	double keff;
	// double conj_keff;
	FLUX *flux;
	//FLUX *conj_flux;
} SSOL;

SSOL *ssol_create(MAPPER *mapper);

void ssol_free(SSOL *ssol);

void ssol_fprintf(const SSOL *ssol, FILE *stream);

#endif
