#ifndef MAPPER_H
#define MAPPER_H

#include<stddef.h>
#include"sconf.h"

typedef struct
{
	size_t xi;
	size_t yi;
	size_t zi;
} XYZ_IDX;

// mapper should contains a checker!!!
typedef struct
{
	size_t xm_size;
	size_t ym_size;
	size_t zm_size;
	size_t rt_size;
	XYZ_IDX *one2three;
	size_t ***three2one; //indexed [k][j][i]
} MAPPER;

MAPPER *mapper_create(const SCONF *sconf);

void mapper_free(MAPPER *mapper);

XYZ_IDX mapper_get3Didx(const MAPPER *mapper, size_t idx1D);

size_t mapper_get1Didx(const MAPPER *mapper, size_t i, size_t j, size_t k);

void mapper_fprintf(const MAPPER *mapper, FILE *stream);

#endif
