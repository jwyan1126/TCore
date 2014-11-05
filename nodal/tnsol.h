#ifndef TNSOL_H
#define TNSOL_H

#include<stddef.h>
#include<stdio.h>

typedef struct
{
	size_t eg_size;
	double keff;
	int bdy;
	double *Dgi;
	double *Dgj;
	double dui;
	double duj;
	double *vsfgi;
	double *vsfgj;
	double *phigi;
	double *phigj;
	double *srgi;
	double *srgj;
	double **ssgi;
	double **ssgj;
	double *chigi;
	double *chigj;
	double *lgi0;
	double *lgj0;
	double *lgi1;
	double *lgj1;
	double *lgi2;
	double *lgj2;
	double *adfgi;
	double *adfgj;
	double *J; // result
	double *agi1;
	double *agj1;
	double *agi2;
	double *agj2;
	double *agi3;
	double *agj3;
	double *agi4;
	double *agj4;
} TNSOL;

TNSOL *tnsol_create(size_t eg_size);

void tnsol_free(TNSOL *tnsol);

void tnsol_coef_fprintf(const TNSOL *tn, FILE *stream);

#endif
