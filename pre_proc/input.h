#ifndef INPUT_H
#define INPUT_H

#include"mtrllib.h"
#include<stddef.h>

typedef struct
{
	size_t eg_size;
	size_t xm_span_size;
	size_t ym_span_size;
	size_t zm_span_size;
	double *xspan_len;
	double *yspan_len;
	double *zspan_len;
	size_t *xspan_subdiv;
	size_t *yspan_subdiv;
	size_t *zspan_subdiv;

	// 0 ref
	// 1 zero-flux
	// 2 vac
	int xl_bdy;
	int xr_bdy;
	int yl_bdy;
	int yr_bdy;
	int zl_bdy;
	int zr_bdy;

	int ***mtrl_set;
	MTRLLIB *mtrllib;

	//transient setting
	size_t pcs_size;
	double *nvel;
	double *lambdas;
	double *betas;
	double tau;
	int steps;
} INPUT;

INPUT *input_create(const char *path);

void input_free(INPUT *input);

void input_fprintf(const INPUT *input, FILE *stream);

#endif
