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
	input->xr_bdy = 2;
	input->yl_bdy = 0;
	input->yr_bdy = 2;
	input->zl_bdy = 0;
	input->zr_bdy = 0;

	input->mtrl_set = calloc(input->xm_span_size, sizeof(int **));
	for(size_t i=0; i<input->xm_span_size; ++i){
		input->mtrl_set[i] = calloc(input->ym_span_size, sizeof(int *));
		for(size_t j=0; j<input->ym_span_size; ++j)
			input->mtrl_set[i][j] = calloc(input->zm_span_size, sizeof(int));
	}

	input->mtrl_set[0][0][0] = 2;
	input->mtrl_set[1][0][0] = 3;
	input->mtrl_set[2][0][0] = 1;
	input->mtrl_set[3][0][0] = 4;
	input->mtrl_set[4][0][0] = 2;
	input->mtrl_set[5][0][0] = 3;
	input->mtrl_set[6][0][0] = 1;
	input->mtrl_set[7][0][0] = 1;
	input->mtrl_set[8][0][0] = 5;

	input->mtrl_set[0][1][0] = 3;
	input->mtrl_set[1][1][0] = 2;
	input->mtrl_set[2][1][0] = 4;
	input->mtrl_set[3][1][0] = 1;
	input->mtrl_set[4][1][0] = 3;
	input->mtrl_set[5][1][0] = 1;
	input->mtrl_set[6][1][0] = 3;
	input->mtrl_set[7][1][0] = 1;
	input->mtrl_set[8][1][0] = 5;

	input->mtrl_set[0][2][0] = 1;
	input->mtrl_set[1][2][0] = 4;
	input->mtrl_set[2][2][0] = 2;
	input->mtrl_set[3][2][0] = 3;
	input->mtrl_set[4][2][0] = 1;
	input->mtrl_set[5][2][0] = 3;
	input->mtrl_set[6][2][0] = 1;
	input->mtrl_set[7][2][0] = 1;
	input->mtrl_set[8][2][0] = 5;

	input->mtrl_set[0][3][0] = 4;
	input->mtrl_set[1][3][0] = 1;
	input->mtrl_set[2][3][0] = 3;
	input->mtrl_set[3][3][0] = 2;
	input->mtrl_set[4][3][0] = 4;
	input->mtrl_set[5][3][0] = 3;
	input->mtrl_set[6][3][0] = 1;
	input->mtrl_set[7][3][0] = 5;
	input->mtrl_set[8][3][0] = 5;

	
	input->mtrl_set[0][4][0] = 2;
	input->mtrl_set[1][4][0] = 3;
	input->mtrl_set[2][4][0] = 1;
	input->mtrl_set[3][4][0] = 4;
	input->mtrl_set[4][4][0] = 2;
	input->mtrl_set[5][4][0] = 1;
	input->mtrl_set[6][4][0] = 1;
	input->mtrl_set[7][4][0] = 5;
	input->mtrl_set[8][4][0] = -1;

	input->mtrl_set[0][5][0] = 3;
	input->mtrl_set[1][5][0] = 1;
	input->mtrl_set[2][5][0] = 3;
	input->mtrl_set[3][5][0] = 3;
	input->mtrl_set[4][5][0] = 1;
	input->mtrl_set[5][5][0] = 1;
	input->mtrl_set[6][5][0] = 5;
	input->mtrl_set[7][5][0] = 5;
	input->mtrl_set[8][5][0] = -1;

	input->mtrl_set[0][6][0] = 1;
	input->mtrl_set[1][6][0] = 3;
	input->mtrl_set[2][6][0] = 1;
	input->mtrl_set[3][6][0] = 1;
	input->mtrl_set[4][6][0] = 1;
	input->mtrl_set[5][6][0] = 5;
	input->mtrl_set[6][6][0] = 5;
	input->mtrl_set[7][6][0] = -1;
	input->mtrl_set[8][6][0] = -1;

	input->mtrl_set[0][7][0] = 1;
	input->mtrl_set[1][7][0] = 1;
	input->mtrl_set[2][7][0] = 1;
	input->mtrl_set[3][7][0] = 5;
	input->mtrl_set[4][7][0] = 5;
	input->mtrl_set[5][7][0] = 5;
	input->mtrl_set[6][7][0] = -1;
	input->mtrl_set[7][7][0] = -1;
	input->mtrl_set[8][7][0] = -1;

	input->mtrl_set[0][8][0] = 5;
	input->mtrl_set[1][8][0] = 5;
	input->mtrl_set[2][8][0] = 5;
	input->mtrl_set[3][8][0] = 5;
	input->mtrl_set[4][8][0] = -1;
	input->mtrl_set[5][8][0] = -1;
	input->mtrl_set[6][8][0] = -1;
	input->mtrl_set[7][8][0] = -1;
	input->mtrl_set[8][8][0] = -1;

	for(size_t i=0; i<xm_span_size; ++i){
		input->xspan_len[i] = 15.0;
		input->xspan_subdiv[i] = 1;
	}
	for(size_t j=0; j<ym_span_size; ++j){
		input->yspan_len[j] = 15.0;
		input->yspan_subdiv[j] = 1;
	}
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
	dcoef[0] = 1.8440; dcoef[1] = 0.4284;
	sa[0] = 0.00607; sa[1] = 0.05946;
	vsf[0] = 0.004556; vsf[1] = 0.072540;
	ss[1][0] = 0.01874;
	adfxr[0] = 0.9623; adfxr[1] = 1.4510;
	adfyr[0] = 0.9623; adfyr[1] = 1.4510;
	adfxl[0] = 0.9623; adfxl[1] = 1.4510;
	adfyl[0] = 0.9623; adfyl[1] = 1.4510;
	MTRL *m1 = mtrl_create(1,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);
	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.8580; dcoef[1] = 0.4283;
	sa[0] = 0.00804; sa[1] = 0.07416;
	vsf[0] = 0.004566; vsf[1] = 0.075580;
	ss[1][0] = 0.01772;
	adfxr[0] = 1.0150; adfxr[1] = 1.8880;
	adfyr[0] = 1.0150; adfyr[1] = 1.8880;
	adfxl[0] = 0.8955; adfxl[1] = 0.6492;
	adfyl[0] = 0.8955; adfyl[1] = 0.6492;
	MTRL *m2 = mtrl_create(2,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);

	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.8440; dcoef[1] = 0.4284;
	sa[0] = 0.00608; sa[1] = 0.05946;
	vsf[0] = 0.003796; vsf[1] = 0.065950;
	ss[1][0] = 0.01874;
	adfxr[0] = 0.9625; adfxr[1] = 1.4510;
	adfyr[0] = 0.9625; adfyr[1] = 1.4510;
	adfxl[0] = 0.9625; adfxl[1] = 1.4510;
	adfyl[0] = 0.9625; adfyl[1] = 1.4510;
	MTRL *m3 = mtrl_create(3,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);

	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.8580; dcoef[1] = 0.4283;
	sa[0] = 0.00804; sa[1] = 0.07415;
	vsf[0] = 0.003804; vsf[1] = 0.068700;
	ss[1][0] = 0.01772;
	adfxr[0] = 1.0160; adfxr[1] = 1.8890;
	adfyr[0] = 0.8949; adfyr[1] = 0.6488;
	adfxl[0] = 0.8949; adfxl[1] = 0.6488;
	adfyl[0] = 1.0160; adfyl[1] = 1.8890;
	MTRL *m4 = mtrl_create(4,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);

	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 2.0000; dcoef[1] = 0.3000;
	sa[0] = 0.0; sa[1] = 0.01;
	vsf[0] = 0.0; vsf[1] = 0.0;
	ss[1][0] = 0.04;
	adfxr[0] = 1.0; adfxr[1] = 1.0;
	adfyr[0] = 1.0; adfyr[1] = 1.0;
	adfxl[0] = 1.0; adfxl[1] = 1.0;
	adfyl[0] = 1.0; adfyl[1] = 1.0;
	MTRL *m5 = mtrl_create(5,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);

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
