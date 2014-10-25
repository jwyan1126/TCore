#ifndef MTRL_H
#define MTRL_H

#include<stddef.h>
#include<stdio.h>

typedef struct
{
	int mtrl_id; // Material ID
	size_t eg_size; // Num of energy groups
	double *chi; // Fission spectral
	double *dcoef; // Diffusion coefficient
	double *sa; // Sigma_a
	double *sr; // Sigma_r
	double *vsf; // vSigma_f
	double **ss; // Sigma_s
	double *adfxl; // Left ADF in x direction
	double *adfxr; // Right ADF in x direction
	double *adfyl; // Left ADF in y direction
	double *adfyr; // Right ADF in y direction
	double *adfzl; // Left ADF in z direction
	double *adfzr; // Right ADF in z direction
} MTRL;

MTRL *mtrl_create(int mtrl_id,
		  size_t eg_size,
		  double *chi,
		  double *dcoef,
		  double *sa,
		  double *vsf,
		  double **ss,
		  double *adfxl,
		  double *adfxr,
		  double *adfyl,
		  double *adfyr,
		  double *adfzl,
		  double *adfzr);

void mtrl_free(MTRL *m);

void mtrl_fprintf(const MTRL *m, FILE *stream);

int mtrl_get_id(const MTRL *m);

double mtrl_get_chi(const MTRL *m, size_t g);

double mtrl_get_dcoef(const MTRL *m, size_t g);

double mtrl_get_sa(const MTRL *m, size_t g);

double mtrl_get_sr(const MTRL *m, size_t g);

double mtrl_get_vsf(const MTRL *m, size_t g);

double mtrl_get_ss(const MTRL *m, size_t g, size_t from_g);

double mtrl_get_adfxl(const MTRL *m, size_t g);

double mtrl_get_adfxr(const MTRL *m, size_t g);

double mtrl_get_adfyl(const MTRL *m, size_t g);

double mtrl_get_adfyr(const MTRL *m, size_t g);

double mtrl_get_adfzl(const MTRL *m, size_t g);

double mtrl_get_adfzr(const MTRL *m, size_t g);

#endif
