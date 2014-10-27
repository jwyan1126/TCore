#include"input.h"
#include<stdlib.h>

INPUT *input_create(const char *path)
{
	INPUT *input = malloc(sizeof(INPUT));
	// test
	size_t eg_size = 2;
	size_t xm_span_size = 4;
	size_t ym_span_size = 4;
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

	input->xl_bdy = 0.0;
	input->xr_bdy = 0.5;
	input->yl_bdy = 0.0;
	input->yr_bdy = 0.5;
	input->zl_bdy = 0.0;
	input->zr_bdy = 0.0;

	input->mtrl_set = calloc(input->xm_span_size, sizeof(int **));
	for(size_t i=0; i<input->xm_span_size; ++i){
		input->mtrl_set[i] = calloc(input->ym_span_size, sizeof(int *));
		for(size_t j=0; j<input->ym_span_size; ++j)
			input->mtrl_set[i][j] = calloc(input->zm_span_size, sizeof(int));
	}
	input->mtrl_set[0][0][0] = 1;
	input->mtrl_set[1][0][0] = 2;
	input->mtrl_set[2][0][0] = 1;
	input->mtrl_set[3][0][0] = 3;
	input->mtrl_set[0][1][0] = 2;
	input->mtrl_set[1][1][0] = 1;
	input->mtrl_set[2][1][0] = 2;
	input->mtrl_set[3][1][0] = 3;
	input->mtrl_set[0][2][0] = 1;
	input->mtrl_set[1][2][0] = 2;
	input->mtrl_set[2][2][0] = 1;
	input->mtrl_set[3][2][0] = 3;
	input->mtrl_set[0][3][0] = 3;
	input->mtrl_set[1][3][0] = 3;
	input->mtrl_set[2][3][0] = 3;
	input->mtrl_set[3][3][0] = 3;

	for(size_t i=0; i<xm_span_size; ++i){
		input->xspan_len[i] = 8.0;
		input->xspan_subdiv[i] = 2;
	}
	input->xspan_len[0] = 4.0;
	input->xspan_subdiv[0] = 1;
	for(size_t j=0; j<ym_span_size; ++j){
		input->yspan_len[j] = 8.0;
		input->yspan_subdiv[j] = 2;
	}
	input->yspan_len[0] = 4.0;
	input->yspan_subdiv[0] = 1;
	for(size_t k=0; k<zm_span_size; ++k){
		input->zspan_len[k] = 1.0;
		input->zspan_subdiv[k] = 1;
	}
	input->zspan_len[0] = 1.0;
	input->zspan_subdiv[0] = 1;

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
	dcoef[0] = 1.3880; dcoef[1] = 0.3842;
	sa[0] = 0.00971; sa[1] = 0.121800;
	vsf[0] = 0.006499; vsf[1] = 0.142900;
	ss[1][0] = 0.01422;
	adfxr[0] = 0.9976; adfxr[1] = 0.7996;
	adfyr[0] = 0.9976; adfyr[1] = 0.7996;
	adfxl[0] = 1.0020; adfxl[1] = 1.1390;
	adfyl[0] = 1.0020; adfyl[1] = 1.1390;
	MTRL *m1 = mtrl_create(1,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);
	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.4470; dcoef[1] = 0.3748;
	sa[0] = 0.009447; sa[1] = 0.089290;
	vsf[0] = 0.006505; vsf[1] = 0.132500;
	ss[1][0] = 0.01730;
	adfxr[0] = 0.9939; adfxr[1] = 1.1180;
	adfyr[0] = 0.9939; adfyr[1] = 1.1180;
	adfxl[0] = 1.0040; adfxl[1] = 0.9222;
	adfyl[0] = 1.0040; adfyl[1] = 0.9222;
	MTRL *m2 = mtrl_create(2,eg_size,chi,dcoef,sa,vsf,ss,adfxl,adfxr,adfyl,adfyr,NULL,NULL);
	chi[0] = 1.0; chi[1] = 0.0;
	dcoef[0] = 1.5450; dcoef[1] = 0.3126;
	sa[0] = 0.000444; sa[1] = 0.008736;
	vsf[0] = 0.0; vsf[1] = 0.0;
	ss[1][0] = 0.02838;
	MTRL *m3 = mtrl_create(3,eg_size,chi,dcoef,sa,vsf,ss,NULL,NULL,NULL,NULL,NULL,NULL);
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
	fprintf(stream, "%4g\t%4g\n", input->xl_bdy, input->xr_bdy);
	fprintf(stream, "YL_BDY\tYR_BDY:\n");
	fprintf(stream, "%4g\t%4g\n", input->yl_bdy, input->yr_bdy);
	fprintf(stream, "ZL_BDY\tZR_BDY:\n");
	fprintf(stream, "%4g\t%4g\n", input->zl_bdy, input->zr_bdy);
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
