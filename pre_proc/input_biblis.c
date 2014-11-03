#include"input.h"
#include<stdlib.h>

INPUT *input_create(const char *path)
{
	INPUT *input = malloc(sizeof(INPUT));
	// test
	size_t eg_size = 2;
	size_t xm_span_size = 9;
	size_t ym_span_size = 9;
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

	input->xl_bdy = 0;
	input->xr_bdy = 1;
	input->yl_bdy = 0;
	input->yr_bdy = 1;
	input->zl_bdy = 0;
	input->zr_bdy = 0;

	input->mtrl_set = calloc(input->xm_span_size, sizeof(int **));
	for(size_t i=0; i<input->xm_span_size; ++i){
		input->mtrl_set[i] = calloc(input->ym_span_size, sizeof(int *));
		for(size_t j=0; j<input->ym_span_size; ++j)
			input->mtrl_set[i][j] = calloc(input->zm_span_size, sizeof(int));
	}

	input->mtrl_set[0][0][0] = 1;
	input->mtrl_set[1][0][0] = 8;
	input->mtrl_set[2][0][0] = 2;
	input->mtrl_set[3][0][0] = 6;
	input->mtrl_set[4][0][0] = 1;
	input->mtrl_set[5][0][0] = 7;
	input->mtrl_set[6][0][0] = 1;
	input->mtrl_set[7][0][0] = 4;
	input->mtrl_set[8][0][0] = 3;

	input->mtrl_set[0][1][0] = 8;
	input->mtrl_set[1][1][0] = 1;
	input->mtrl_set[2][1][0] = 8;
	input->mtrl_set[3][1][0] = 2;
	input->mtrl_set[4][1][0] = 8;
	input->mtrl_set[5][1][0] = 1;
	input->mtrl_set[6][1][0] = 1;
	input->mtrl_set[7][1][0] = 4;
	input->mtrl_set[8][1][0] = 3;

	input->mtrl_set[0][2][0] = 2;
	input->mtrl_set[1][2][0] = 8;
	input->mtrl_set[2][2][0] = 1;
	input->mtrl_set[3][2][0] = 8;
	input->mtrl_set[4][2][0] = 2;
	input->mtrl_set[5][2][0] = 7;
	input->mtrl_set[6][2][0] = 1;
	input->mtrl_set[7][2][0] = 4;
	input->mtrl_set[8][2][0] = 3;

	input->mtrl_set[0][3][0] = 6;
	input->mtrl_set[1][3][0] = 2;
	input->mtrl_set[2][3][0] = 8;
	input->mtrl_set[3][3][0] = 2;
	input->mtrl_set[4][3][0] = 8;
	input->mtrl_set[5][3][0] = 1;
	input->mtrl_set[6][3][0] = 8;
	input->mtrl_set[7][3][0] = 4;
	input->mtrl_set[8][3][0] = 3;

	
	input->mtrl_set[0][4][0] = 1;
	input->mtrl_set[1][4][0] = 8;
	input->mtrl_set[2][4][0] = 2;
	input->mtrl_set[3][4][0] = 8;
	input->mtrl_set[4][4][0] = 2;
	input->mtrl_set[5][4][0] = 5;
	input->mtrl_set[6][4][0] = 4;
	input->mtrl_set[7][4][0] = 3;
	input->mtrl_set[8][4][0] = 3;

	input->mtrl_set[0][5][0] = 7;
	input->mtrl_set[1][5][0] = 1;
	input->mtrl_set[2][5][0] = 7;
	input->mtrl_set[3][5][0] = 1;
	input->mtrl_set[4][5][0] = 5;
	input->mtrl_set[5][5][0] = 4;
	input->mtrl_set[6][5][0] = 4;
	input->mtrl_set[7][5][0] = 3;
	input->mtrl_set[8][5][0] = -1;

	input->mtrl_set[0][6][0] = 1;
	input->mtrl_set[1][6][0] = 1;
	input->mtrl_set[2][6][0] = 1;
	input->mtrl_set[3][6][0] = 8;
	input->mtrl_set[4][6][0] = 4;
	input->mtrl_set[5][6][0] = 4;
	input->mtrl_set[6][6][0] = 3;
	input->mtrl_set[7][6][0] = 3;
	input->mtrl_set[8][6][0] = -1;

	input->mtrl_set[0][7][0] = 4;
	input->mtrl_set[1][7][0] = 4;
	input->mtrl_set[2][7][0] = 4;
	input->mtrl_set[3][7][0] = 4;
	input->mtrl_set[4][7][0] = 3;
	input->mtrl_set[5][7][0] = 3;
	input->mtrl_set[6][7][0] = 3;
	input->mtrl_set[7][7][0] = -1;
	input->mtrl_set[8][7][0] = -1;

	input->mtrl_set[0][8][0] = 3;
	input->mtrl_set[1][8][0] = 3;
	input->mtrl_set[2][8][0] = 3;
	input->mtrl_set[3][8][0] = 3;
	input->mtrl_set[4][8][0] = 3;
	input->mtrl_set[5][8][0] = -1;
	input->mtrl_set[6][8][0] = -1;
	input->mtrl_set[7][8][0] = -1;
	input->mtrl_set[8][8][0] = -1;

	for(size_t i=0; i<xm_span_size; ++i){
		input->xspan_len[i] = 23.1226;
		input->xspan_subdiv[i] = 2;
	}
	input->xspan_len[0] = 23.1226 / 2.0;
	input->xspan_subdiv[0] = 1;
	for(size_t j=0; j<ym_span_size; ++j){
		input->yspan_len[j] = 23.1226;
		input->yspan_subdiv[j] = 2;
	}
	input->yspan_len[0] = 23.1226 / 2.0;
	input->yspan_subdiv[0] = 1;
	for(size_t k=0; k<zm_span_size; ++k){
		input->zspan_len[k] = 1.0;
		input->zspan_subdiv[k] = 1;
	}

	input->mtrllib = mtrllib_create();
	double *chi = calloc(2,sizeof(double));
	double *dcoef = calloc(2,sizeof(double));
	double *sa = calloc(2,sizeof(double));
	double *vsf = calloc(2,sizeof(double));
	double **ss = calloc(2,sizeof(double *));
	double *adfxl = calloc(2,sizeof(double));
	double *adfxr = calloc(2,sizeof(double));
	double *adfyl = calloc(2,sizeof(double));
	double *adfyr = calloc(2,sizeof(double));
	double *adfzl = calloc(2,sizeof(double));
	double *adfzr = calloc(2,sizeof(double));
	for(size_t i=0; i<2; ++i)
		ss[i] = calloc(2,sizeof(double));

	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.436e0; dcoef[1] = 3.6350e-1;
	sa[0] = 9.5042e-3; sa[1] = 7.5058e-2;
	vsf[0] = 5.8708e-3; vsf[1] = 9.6067e-2;
	ss[1][0] = 1.7754e-2;
	adfxr[0] = 1.0; adfxr[1] = 1.0;
	adfyr[0] = 1.0; adfyr[1] = 1.0;
	adfxl[0] = 1.0; adfxl[1] = 1.0;
	adfyl[0] = 1.0; adfyl[1] = 1.0;
	MTRL *m1 = mtrl_create(1,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);
	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.4366e0; dcoef[1] = 3.6360e-1;
	sa[0] = 9.6785e-3; sa[1] = 7.8436e-2;
	vsf[0] = 6.1908e-3; vsf[1] = 1.0358e-1;
	ss[1][0] = 1.7621e-2;
	adfxr[0] = 1.0; adfxr[1] = 1.0;
	adfyr[0] = 1.0; adfyr[1] = 1.0;
	adfxl[0] = 1.0; adfxl[1] = 1.0;
	adfyl[0] = 1.0; adfyl[1] = 1.0;
	MTRL *m2 = mtrl_create(2,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);

	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.32e0; dcoef[1] = 2.772e-1;
	sa[0] = 2.6562e-3; sa[1] = 7.1596e-2;
	vsf[0] = 0.0; vsf[1] = 0.0;
	ss[1][0] = 2.3106e-2;
	adfxr[0] = 1.0; adfxr[1] = 1.0;
	adfyr[0] = 1.0; adfyr[1] = 1.0;
	adfxl[0] = 1.0; adfxl[1] = 1.0;
	adfyl[0] = 1.0; adfyl[1] = 1.0;
	MTRL *m3 = mtrl_create(3,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);

	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.4389e0; dcoef[1] = 3.638e-1;
	sa[0] = 1.0363e-2; sa[1] = 9.1408e-2;
	vsf[0] = 7.4527e-3; vsf[1] = 1.3236e-1;
	ss[1][0] = 1.71010e-2;
	adfxr[0] = 1.0; adfxr[1] = 1.0;
	adfyr[0] = 1.0; adfyr[1] = 1.0;
	adfxl[0] = 1.0; adfxl[1] = 1.0;
	adfyl[0] = 1.0; adfyl[1] = 1.0;
	MTRL *m4 = mtrl_create(4,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);

	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.4381e0; dcoef[1] = 3.665e-1;
	sa[0] = 1.0003e-2; sa[1] = 8.4828e-2;
	vsf[0] = 6.1908e-3; vsf[1] = 1.0358e-1;
	ss[1][0] = 1.7290e-2;
	adfxr[0] = 1.0; adfxr[1] = 1.0;
	adfyr[0] = 1.0; adfyr[1] = 1.0;
	adfxl[0] = 1.0; adfxl[1] = 1.0;
	adfyl[0] = 1.0; adfyl[1] = 1.0;
	MTRL *m5 = mtrl_create(5,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);

	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.4385e0; dcoef[1] = 3.665e-1;
	sa[0] = 1.0132e-2; sa[1] = 8.7314e-2;
	vsf[0] = 6.4285e-3; vsf[1] = 1.0911e-1;
	ss[1][0] = 1.7192e-2;
	adfxr[0] = 1.0; adfxr[1] = 1.0;
	adfyr[0] = 1.0; adfyr[1] = 1.0;
	adfxl[0] = 1.0; adfxl[1] = 1.0;
	adfyl[0] = 1.0; adfyl[1] = 1.0;
	MTRL *m6 = mtrl_create(6,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);

	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.4389e0; dcoef[1] = 3.679e-1;
	sa[0] = 1.0165e-2; sa[1] = 8.8024e-2;
	vsf[0] = 6.1908e-3; vsf[1] = 1.0358e-1;
	ss[1][0] = 1.7125e-2;
	adfxr[0] = 1.0; adfxr[1] = 1.0;
	adfyr[0] = 1.0; adfyr[1] = 1.0;
	adfxl[0] = 1.0; adfxl[1] = 1.0;
	adfyl[0] = 1.0; adfyl[1] = 1.0;
	MTRL *m7 = mtrl_create(7,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);

	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.4393e0; dcoef[1] = 3.68e-1;
	sa[0] = 1.0294e-2; sa[1] = 9.051e-2;
	vsf[0] = 6.4285e-3; vsf[1] = 1.0911e-1;
	ss[1][0] = 1.7027e-2;
	adfxr[0] = 1.0; adfxr[1] = 1.0;
	adfyr[0] = 1.0; adfyr[1] = 1.0;
	adfxl[0] = 1.0; adfxl[1] = 1.0;
	adfyl[0] = 1.0; adfyl[1] = 1.0;
	MTRL *m8 = mtrl_create(8,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);
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
	mtrllib_add(input->mtrllib, m3);
	mtrllib_add(input->mtrllib, m4);
	mtrllib_add(input->mtrllib, m5);
	mtrllib_add(input->mtrllib, m6);
	mtrllib_add(input->mtrllib, m7);
	mtrllib_add(input->mtrllib, m8);
	
	// ...
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
