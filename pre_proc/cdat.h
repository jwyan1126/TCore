#ifndef CDAT_H
#define CDAT_H

#include"sconf.h"
#include<stddef.h>
#include<stdio.h>

typedef struct
{
	size_t xsize;
	size_t ysize;
	size_t zsize;
	int ***checker;
	double ***data;
} CDAT3;

CDAT3 *cdat3_create(SCONF *sconf);

void cdat3_free(CDAT3 *dat);

double cdat3_get_val(const CDAT3 *dat, size_t i, size_t j, size_t k);

void cdat3_set_val(CDAT3 *dat, size_t i, size_t j, size_t k, double val);

void cdat3_copy(CDAT3 *tar_dat, const CDAT3 *src_dat);

typedef struct
{
	size_t xsize;
	size_t ysize;
	size_t zsize;
	size_t gsize;
	int ***checker;
	double ****data;
} CDAT4;

CDAT4 *cdat4_create(SCONF *sconf);

void cdat4_free(CDAT4 *dat);

double cdat4_get_val(const CDAT4 *dat, size_t g, size_t i, size_t j, size_t k);

void cdat4_set_val(CDAT4 *dat, size_t g, size_t i, size_t j, size_t k, double val);

void cdat4_copy(CDAT4 *tar_dat, const CDAT4 *src_dat);

typedef struct
{
	size_t xsize;
	size_t ysize;
	size_t zsize;
	size_t gsize;
	int ***checker;
	double *****data;
} CDAT5;

CDAT5 *cdat5_create(SCONF *sconf);

void cdat5_free(CDAT5 *dat);

double cdat5_get_val(const CDAT5 *dat, size_t g, size_t from_g, size_t i, size_t j, size_t k);

void cdat5_set_val(CDAT5 *dat, size_t g, size_t from_g, size_t i, size_t j, size_t k, double val);

void cdat5_copy(CDAT5 *tar_dat, const CDAT5 *src_dat);

#endif
