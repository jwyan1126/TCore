#include"input.h"
#include<stdlib.h>

INPUT *input_create(const char *path)
{
	INPUT *input = malloc(sizeof(INPUT));
	// test
	size_t eg_size = 2;
	size_t xm_span_size = 3;
	size_t ym_span_size = 1;
	size_t zm_span_size = 1;
	input->eg_size = eg_size;
	input->xm_span_size = xm_span_size;
	input->ym_span_size = ym_span_size;
	input->zm_span_size = zm_span_size;
	input->xspan_len = calloc(xm_span_size, sizeof(double));
	input->yspan_len = calloc(ym_span_size, sizeof(double));
	input->zspan_len = calloc(zm_span_size, sizeof(double));
	input->xspan_subdiv = calloc(xm_span_size, sizeof(size_t));
	input->yspan_subdiv = calloc(ym_span_size, sizeof(size_t));
	input->zspan_subdiv = calloc(zm_span_size, sizeof(size_t));

	input->xl_bdy = 1;
	input->xr_bdy = 1;
	input->yl_bdy = 0;
	input->yr_bdy = 0;
	input->zl_bdy = 0;
	input->zr_bdy = 0;

	input->mtrl_set = calloc(input->xm_span_size, sizeof(int **));
	for(size_t i=0; i<input->xm_span_size; ++i){
		input->mtrl_set[i] = calloc(input->ym_span_size, sizeof(int *));
		for(size_t j=0; j<input->ym_span_size; ++j)
			input->mtrl_set[i][j] = calloc(input->zm_span_size, sizeof(int));
	}

	input->mtrl_set[0][0][0] = 1;
	input->mtrl_set[1][0][0] = 2;
	input->mtrl_set[2][0][0] = 1;

	input->xspan_len[0] = 40.0;
	input->xspan_len[1] = 160.0;
	input->xspan_len[2] = 40.0;
	input->xspan_subdiv[0] = 1;
	input->xspan_subdiv[1] = 4;
	input->xspan_subdiv[2] = 1;

	input->yspan_len[0] = 1.0;
	input->yspan_subdiv[0] = 1;

	input->zspan_len[0] = 1.0;
	input->zspan_subdiv[0] = 1;

	input->mtrllib = mtrllib_create();
	double *chi = calloc(eg_size,sizeof(double));
	double *dcoef = calloc(eg_size,sizeof(double));
	double *sa = calloc(eg_size,sizeof(double));
	double *vsf = calloc(eg_size,sizeof(double));
	double **ss = calloc(eg_size,sizeof(double *));
	double *adfxl = calloc(eg_size,sizeof(double));
	double *adfxr = calloc(eg_size,sizeof(double));
	double *adfyl = calloc(eg_size,sizeof(double));
	double *adfyr = calloc(eg_size,sizeof(double));
	double *adfzl = calloc(eg_size,sizeof(double));
	double *adfzr = calloc(eg_size,sizeof(double));
	for(size_t i=0; i<eg_size; ++i)
		ss[i] = calloc(eg_size,sizeof(double));

	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.5; dcoef[1] = 0.5;
	sa[0] = 0.026-0.015; sa[1] = 0.18;
	vsf[0] = 0.01; vsf[1] = 0.2;
	ss[1][0] = 0.015;
	adfxr[0] = 1.0; adfxr[1] = 1.0;
	adfyr[0] = 1.0; adfyr[1] = 1.0;
	adfxl[0] = 1.0; adfxl[1] = 1.0;
	adfyl[0] = 1.0; adfyl[1] = 1.0;
	MTRL *m1 = mtrl_create(1,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);

	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.0; dcoef[1] = 0.5;
	sa[0] = 0.02-0.01; sa[1] = 0.08;
	vsf[0] = 0.005; vsf[1] = 0.099;
	ss[1][0] = 0.01;
	adfxr[0] = 1.0; adfxr[1] = 1.0;
	adfyr[0] = 1.0; adfyr[1] = 1.0;
	adfxl[0] = 1.0; adfxl[1] = 1.0;
	adfyl[0] = 1.0; adfyl[1] = 1.0;
	MTRL *m2 = mtrl_create(2,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);
	

	free(adfxl);
	free(adfxr);
	free(adfyl);
	free(adfyr);
	free(adfzl);
	free(adfzr);
	free(chi);
	free(dcoef);
	free(sa);
	free(vsf);
	for(size_t i=0; i<2; ++i)
		free(ss[i]);
	free(ss);
	mtrllib_add(input->mtrllib, m1);
	mtrllib_add(input->mtrllib, m2);
	
	
	size_t pcs_size = 6;
	input->pcs_size = pcs_size;
	input->nvel = malloc(eg_size * sizeof(double));
	input->lambdas = malloc(pcs_size * sizeof(double));
	input->betas = malloc(pcs_size * sizeof(double));
	input->nvel[0] = 1.0e7;
	input->nvel[1] = 3.0e5;
	input->lambdas[0] = 0.0124;
	input->lambdas[1] = 0.0305;
	input->lambdas[2] = 0.1110;
	input->lambdas[3] = 0.3010;
	input->lambdas[4] = 1.1400;
	input->lambdas[5] = 3.0100;
	input->betas[0] = 0.00025;
	input->betas[1] = 0.00164;
	input->betas[2] = 0.00147;
	input->betas[3] = 0.00296;
	input->betas[4] = 0.00086;
	input->betas[5] = 0.00032;
	input->tau = 0.1;
	input->steps = 100;
	return input;
}

void input_free(INPUT *input)
{
	free(input->xspan_len);
	free(input->yspan_len);
	free(input->zspan_len);
	free(input->xspan_subdiv);
	free(input->yspan_subdiv);
	free(input->zspan_subdiv);
	for(size_t i=0; i<input->xm_span_size; ++i){
		for(size_t j=0; j<input->ym_span_size; ++j)
			free(input->mtrl_set[i][j]);
		free(input->mtrl_set[i]);
	}
	free(input->mtrl_set);
	free(input);
}

void input_fprintf(const INPUT *input, FILE *stream)
{
	fprintf(stream, "EG_SIZE:\n");
	fprintf(stream, "%4zd\n", input->eg_size);
	fprintf(stream, "XM_SPAN_SIZE\tYM_SPAN_SIZE\tZM_SPAN_SIZE:\n");
	fprintf(stream, "%4zd\t%4zd\t%4zd\n", input->xm_span_size, input->ym_span_size, input->zm_span_size);
	fprintf(stream, "XSPAN_LEN:\n");
	for(size_t i=0; i<input->xm_span_size; ++i)
		fprintf(stream, "%4g\t", input->xspan_len[i]);
	fprintf(stream, "\n");
	fprintf(stream, "YSPAN_LEN:\n");
	for(size_t i=0; i<input->ym_span_size; ++i)
		fprintf(stream, "%4g\t", input->yspan_len[i]);
	fprintf(stream, "\n");
	fprintf(stream, "ZSPAN_LEN:\n");
	for(size_t i=0; i<input->zm_span_size; ++i)
		fprintf(stream, "%4g\t", input->zspan_len[i]);
	fprintf(stream, "\n");
	fprintf(stream, "XSPAN_SUBDIV:\n");
	for(size_t i=0; i<input->xm_span_size; ++i)
		fprintf(stream, "%4zd\t", input->xspan_subdiv[i]);
	fprintf(stream, "\n");
	fprintf(stream, "YSPAN_SUBDIV:\n");
	for(size_t i=0; i<input->ym_span_size; ++i)
		fprintf(stream, "%4zd\t", input->yspan_subdiv[i]);
	fprintf(stream, "\n");
	fprintf(stream, "ZSPAN_SUBDIV:\n");
	for(size_t i=0; i<input->zm_span_size; ++i)
		fprintf(stream, "%4zd\t", input->zspan_subdiv[i]);
	fprintf(stream, "\n");
	fprintf(stream, "XL_BDY\tXR_BDY:\n");
	fprintf(stream, "%4d\t%4d\n", input->xl_bdy, input->xr_bdy);
	fprintf(stream, "YL_BDY\tYR_BDY:\n");
	fprintf(stream, "%4d\t%4d\n", input->yl_bdy, input->yr_bdy);
	fprintf(stream, "ZL_BDY\tZR_BDY:\n");
	fprintf(stream, "%4d\t%4d\n", input->zl_bdy, input->zr_bdy);
	fprintf(stream, "\n");
	fprintf(stream, "MATERIAL SETTING:\n");
	for(size_t k=0; k<input->zm_span_size; ++k){
		fprintf(stream, "z = %zd\n", k);
		for(size_t i=0; i<input->xm_span_size; ++i){
			for(size_t j=0; j<input->ym_span_size; ++j)
				fprintf(stream,"%4d\t", input->mtrl_set[i][j][k]);
			fprintf(stream, "\n");
		}
	}
	fprintf(stream, "\n");
	mtrllib_fprintf(input->mtrllib, stream);
}
