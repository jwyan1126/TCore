#ifndef POWER_H
#define POWER_H

#include"pre_proc/mapper.h"
#include"flux.h"
#include"pre_proc/cdat.h"

typedef struct
{
	size_t xm_size;
	size_t ym_size;
	size_t zm_size;
	size_t rt_size;
	double *data;
	MAPPER *mapper;
} POWER;

POWER *power_create(MAPPER *mapper);

void power_free(POWER *power);

void power_fprintf(const POWER *power, FILE *stream);

void power_normalize(POWER *power);

double power_get_val(const POWER *power, size_t i, size_t j, size_t k);

void power_cal(POWER *power, const FLUX *flux, const CDAT4 *vsf);

double power_sumup(const POWER *power, const CDAT3 *dx, const CDAT3 *dy, const CDAT3 *dz);

#endif
