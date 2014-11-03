#ifndef TNSOL_H
#define TNSOL_H

#include<stddef.h>

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
} TNSOL;

TNSOL *tnsol_create(size_t eg_size);

void tnsol_free(TNSOL *tnsol);

#endif
