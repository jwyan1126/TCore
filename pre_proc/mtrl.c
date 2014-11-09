#include"mtrl.h"
#include<stdio.h>
#include<stdlib.h>

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
		  double *adfzr)
{
	#ifdef DEBUG
	if(mtrl_id < 0){
		fprintf(stderr, "Material ID must be positive.\n");
		exit(-1);
	}
	#endif
	MTRL *mtrl = malloc(sizeof(MTRL));
	mtrl->mtrl_id = mtrl_id;
	mtrl->eg_size = eg_size;
	mtrl->chi = malloc(eg_size * sizeof(double));
	mtrl->dcoef = malloc(eg_size * sizeof(double));
	mtrl->sa = malloc(eg_size * sizeof(double));
	mtrl->sr = malloc(eg_size * sizeof(double));
	mtrl->vsf = malloc(eg_size * sizeof(double));
	mtrl->ss = malloc(eg_size * sizeof(double *));
	mtrl->adfxl = malloc(eg_size * sizeof(double));
	mtrl->adfxr = malloc(eg_size * sizeof(double));
	mtrl->adfyl = malloc(eg_size * sizeof(double));
	mtrl->adfyr = malloc(eg_size * sizeof(double));
	mtrl->adfzl = malloc(eg_size * sizeof(double));
	mtrl->adfzr = malloc(eg_size * sizeof(double));
	for(size_t g=0; g<eg_size; ++g)
		mtrl->ss[g] = malloc(eg_size * sizeof(double));
	for(size_t g=0; g<eg_size; ++g){
		mtrl->chi[g] = chi[g];
		mtrl->dcoef[g] = dcoef[g];
		mtrl->sa[g] = sa[g];
		for(size_t bg=0; bg<eg_size; ++bg)
			mtrl->ss[g][bg] = ss[g][bg];
		mtrl->vsf[g] = vsf[g];
		mtrl->adfxl[g] = (adfxl == NULL) ? 1.0 : adfxl[g];
		mtrl->adfxr[g] = (adfxr == NULL) ? 1.0 : adfxr[g];
		mtrl->adfyl[g] = (adfyl == NULL) ? 1.0 : adfyl[g];
		mtrl->adfyr[g] = (adfyr == NULL) ? 1.0 : adfyr[g];
		mtrl->adfzl[g] = (adfzl == NULL) ? 1.0 : adfzl[g];
		mtrl->adfzr[g] = (adfzr == NULL) ? 1.0 : adfzr[g];
		mtrl_update_sr(mtrl);
	}
	return mtrl;
}

void mtrl_free(MTRL *m)
{
	free(m->chi);
	free(m->dcoef);
	free(m->sa);
	free(m->sr);
	free(m->vsf);
	size_t eg_size = m->eg_size;
	for(size_t i=0; i<eg_size; ++i)
		free(m->ss[i]);
	free(m->ss);
	free(m);
}

void mtrl_update_sr(MTRL *m)
{
	size_t eg_size = m->eg_size;
	for(size_t g=0; g<eg_size; ++g){
		m->sr[g] = m->sa[g];
		for(size_t bg=0; bg<eg_size; ++bg){
			if(g != bg)
				m->sr[g] += m->ss[bg][g];
		}
	}
}

void mtrl_fprintf(const MTRL *m, FILE *stream)
{
	size_t eg_size = m->eg_size;
	fprintf(stream, "MTRL_ID:%4d\n", m->mtrl_id);
	fprintf(stream, "CHI:\n");
	for(size_t g=0; g<eg_size; ++g)
		fprintf(stream, "%4g\t", m->chi[g]);
	fprintf(stream, "\n");
	
	fprintf(stream, "DCOEF:\n");
	for(size_t g=0; g<eg_size; ++g)
		fprintf(stream, "%4g\t", m->dcoef[g]);
	fprintf(stream, "\n");

	fprintf(stream, "SA:\n");
	for(size_t g=0; g<eg_size; ++g)
		fprintf(stream, "%4g\t", m->sa[g]);
	fprintf(stream, "\n");

	fprintf(stream, "VSF:\n");
	for(size_t g=0; g<eg_size; ++g)
		fprintf(stream, "%4g\t", m->vsf[g]);
	fprintf(stream, "\n");

	fprintf(stream, "SS:\n");
	for(size_t g=0; g<eg_size; ++g){
		for(size_t from_g=0; from_g<eg_size; ++from_g)
			fprintf(stream, "%4g\t", m->ss[g][from_g]);
		fprintf(stream, "\n");
	}
	fprintf(stream, "ADFs:\n");
	fprintf(stream, "X-\tX+\tY-\tY+\tZ-\tZ+\n");
	for(size_t g=0; g<eg_size; ++g)
		fprintf(stream, "%4g\t%4g\t%4g\t%4g\t%4g\t%4g\n",
			m->adfxl[g], m->adfxr[g], m->adfyl[g], m->adfyr[g], m->adfzl[g], m->adfxr[g]);
	fprintf(stream, "\n");
}

inline int mtrl_get_id(const MTRL *m)
{
	return m->mtrl_id;
}

inline double mtrl_get_chi(const MTRL *m, size_t g)
{
	#ifdef DEBUG
	size_t eg_size = m->eg_size;
	if(g >= eg_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->chi[g];
}

inline double mtrl_get_dcoef(const MTRL *m, size_t g)
{
	#ifdef DEBUG
	size_t eg_size = m->eg_size;
	if(g >= eg_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->dcoef[g];
}

inline double mtrl_get_sa(const MTRL *m, size_t g)
{
	
	#ifdef DEBUG
	size_t eg_size = m->eg_size;
	if(g >= eg_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->sa[g];
}

inline double mtrl_get_sr(const MTRL *m, size_t g)
{
	#ifdef DEBUG
	size_t eg_size = m->eg_size;
	if(g >= eg_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->sr[g];
}

inline double mtrl_get_vsf(const MTRL *m, size_t g)
{
	#ifdef DEBUG
	size_t eg_size = m->eg_size;
	if(g >= eg_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->vsf[g];
}

inline double mtrl_get_ss(const MTRL *m, size_t g, size_t from_g)
{
	#ifdef DEBUG
	size_t eg_size = m->eg_size;
	if(g >= eg_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->ss[g][from_g];
}

inline double mtrl_get_adfxl(const MTRL *m, size_t g)
{
	#ifdef DEBUG
	size_t eg_size = m->eg_size;
	if(g >= eg_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->adfxl[g];
}

inline double mtrl_get_adfxr(const MTRL *m, size_t g)
{
	#ifdef DEBUG
	size_t eg_size = m->eg_size;
	if(g >= eg_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->adfxr[g];
}

inline double mtrl_get_adfyl(const MTRL *m, size_t g)
{
	#ifdef DEBUG
	size_t eg_size = m->eg_size;
	if(g >= eg_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->adfyl[g];
}

inline double mtrl_get_adfyr(const MTRL *m, size_t g)
{
	#ifdef DEBUG
	size_t eg_size = m->eg_size;
	if(g >= eg_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->adfyr[g];
}

inline double mtrl_get_adfzl(const MTRL *m, size_t g)
{
	#ifdef DEBUG
	size_t eg_size = m->eg_size;
	if(g >= eg_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->adfzl[g];
}

inline double mtrl_get_adfzr(const MTRL *m, size_t g)
{
	#ifdef DEBUG
	size_t eg_size = m->eg_size;
	if(g >= eg_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->adfzr[g];
}
