#ifndef FLUX_H
#define FLUX_H

#include"pre_proc/mapper.h"

typedef struct
{
	size_t eg_size;
	size_t xm_size;
	size_t ym_size;
	size_t zm_size;
	size_t rt_size;
	double *data;
	MAPPER *mapper;
} FLUX;

FLUX *flux_create(MAPPER *mapper);

void flux_free(FLUX *flux);

void flux_fprintf(const FLUX *flux, FILE *stream);

void flux_normalize(FLUX *flux);

#endif
