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

typedef struct
{
	size_t eg_size;
	size_t xm_size;
	size_t ym_size;
	size_t zm_size;
	size_t rt_size;
	XYZ_IDX *one2three;
	size_t ***three2one; 	// indexed [k][j][i]
	int ***cchecker; 	// indexed [k][j][i]
				// cchecker:
				// NO FILL		00000001
				// FILL BUT NOT BDY	00000010
				// X-			00000100
				// X+			00001000
				// Y-			00010000
				// Y+			00100000
				// Z-			01000000
				// Z+			10000000
	int ***xchecker; 	// indexed [k][j][i]
	int ***ychecker; 	// indexed [i][k][j]
	int ***zchecker; 	// indexed [j][i][k]
				// FILL  0001
				// LFILL 0011
				// RFILL 0101
} MAPPER;

MAPPER *mapper_create(SCONF *sconf);

void mapper_free(MAPPER *mapper);

XYZ_IDX mapper_get3Didx(const MAPPER *mapper, size_t idx1D);

size_t mapper_get1Didx(const MAPPER *mapper, size_t i, size_t j, size_t k);

void mapper_fprintf(const MAPPER *mapper, FILE *stream);

#endif
