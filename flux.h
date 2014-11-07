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

void flux_copy(FLUX *tar_flux, const FLUX *src_flux);

void flux_fprintf(const FLUX *flux, FILE *stream);

void flux_normalize(FLUX *flux);

double flux_sumup(const FLUX *flux);

double flux_get_val(FLUX *flux, size_t g, size_t i, size_t j, size_t k);

#endif
