#ifndef EDAT4_H
#define EDAT4_H

#include<stddef.h>
#include<stdio.h>
#include"sconf.h"
#include"mapper.h"

typedef struct
{
	size_t xsize;
	size_t ysize;
	size_t zsize;
	size_t gsize;
	// FILL  0001
	// LFILL 0011
	// RFILL 0101
	int ***xchecker;
	double ****xdata; // indexed by [z][y][x][g]
	int ***ychecker;
	double ****ydata; // indexed by [x][z][y][g]
	int ***zchecker;
	double ****zdata; // indexed by [y][x][z][g]
} EDAT4;

EDAT4 *edat4_create(MAPPER *mapper);

void edat4_free(EDAT4 *dat);

double edat4_get_xlval(const EDAT4 *dat, size_t g, size_t i, size_t j, size_t k);
double edat4_get_xrval(const EDAT4 *dat, size_t g, size_t i, size_t j, size_t k);
double edat4_get_ylval(const EDAT4 *dat, size_t g, size_t i, size_t j, size_t k);
double edat4_get_yrval(const EDAT4 *dat, size_t g, size_t i, size_t j, size_t k);
double edat4_get_zlval(const EDAT4 *dat, size_t g, size_t i, size_t j, size_t k);
double edat4_get_zrval(const EDAT4 *dat, size_t g, size_t i, size_t j, size_t k);

void edat4_set_xlval(EDAT4 *dat, size_t g, size_t i, size_t j, size_t k, double val);
void edat4_set_xrval(EDAT4 *dat, size_t g, size_t i, size_t j, size_t k, double val);
void edat4_set_ylval(EDAT4 *dat, size_t g, size_t i, size_t j, size_t k, double val);
void edat4_set_yrval(EDAT4 *dat, size_t g, size_t i, size_t j, size_t k, double val);
void edat4_set_zlval(EDAT4 *dat, size_t g, size_t i, size_t j, size_t k, double val);
void edat4_set_zrval(EDAT4 *dat, size_t g, size_t i, size_t j, size_t k, double val);

void edat4_copy(EDAT4 *tar_dat, const EDAT4 *src_dat);

void edat4_fprintf(EDAT4 *dat, size_t g, size_t i, size_t j, size_t k, FILE *stream);

void edat4_set_rand(EDAT4 *dat);

#endif
