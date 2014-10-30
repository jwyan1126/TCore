#ifndef TNINP_H
#define TNINP_H

typedef struct
{
	size_t eg_size;
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
} TNINP;

TNINP *tninp_create(size_t eg_size);

void tninp_free(TNINP *tninp);

#endif
