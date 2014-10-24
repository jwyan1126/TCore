#ifndef MTRL_H
#define MTRL_H

#include<stddef.h>
#include<stdio.h>

//extern size_t EG_SIZE;

typedef struct
{
	int mtrl_id;
	size_t eg_size;
	double *chi;
	double *dcoef;
	double *sa;
	double *sr;
	double *vsf;
	double **ss;
} MTRL;

MTRL *mtrl_create(int mtrl_id,
		  size_t eg_size,
		  double *chi,
		  double *dcoef,
		  double *sa,
		  double *vsf,
		  double **ss);

void mtrl_free(MTRL *m);

void mtrl_fprintf(const MTRL *m, FILE *stream);

int mtrl_get_id(const MTRL *m);

double mtrl_get_chi(const MTRL *m, size_t g);

double mtrl_get_dcoef(const MTRL *m, size_t g);

double mtrl_get_sa(const MTRL *m, size_t g);

double mtrl_get_sr(const MTRL *m, size_t g);

double mtrl_get_vsf(const MTRL *m, size_t g);

double mtrl_get_ss(const MTRL *m, size_t g, size_t from_g);

#endif
