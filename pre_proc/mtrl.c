#include"mtrl.h"
#include<stdio.h>
#include<stdlib.h>

MTRL *mtrl_create(int mtrl_id,
		  double *chi,
		  double *dcoef,
		  double *sa,
		  double *vsf,
		  double **ss)
{
	#ifdef DEBUG
	if(mtrl_id < 0){
		fprintf(stderr, "Material ID must be positive.\n");
		exit(-1);
	}
	#endif
	MTRL *mtrl = malloc(sizeof(MTRL));
	mtrl->mtrl_id = mtrl_id;
	mtrl->chi = malloc(EGSIZE * sizeof(double));
	mtrl->dcoef = malloc(EGSIZE * sizeof(double));
	mtrl->sa = malloc(EGSIZE * sizeof(double));
	mtrl->sr = malloc(EGSIZE * sizeof(double));
	mtrl->vsf = malloc(EGSIZE * sizeof(double));
	mtrl->ss = malloc(EGSIZE * sizeof(double *));
	for(size_t g=0; g<EGSIZE; ++g)
		mtrl->ss[g] = malloc(EGSIZE * sizeof(double));
	for(size_t g=0; g<EGSIZE; ++g){
		mtrl->chi[g] = chi[g];
		mtrl->dcoef[g] = dcoef[g];
		mtrl->sa[g] = sa[g];
		mtrl->sr[g] = sa[g];
		mtrl->vsf[g] = vsf[g];
		for(size_t bg=0; bg<EGSIZE; ++bg){
			mtrl->ss[g][bg] = ss[g][bg];
			if(g != bg)
				mtrl->sr[g] += ss[bg][g];
		}
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
	for(size_t i=0; i<EGSIZE; ++i)
		free(m->ss[i]);
	free(m->ss);
	free(m);
}

inline int mtrl_get_id(const MTRL *m)
{
	return m->mtrl_id;
}

inline double mtrl_get_chi(const MTRL *m, size_t g)
{
	#ifdef DEBUG
	if(g >= EGSIZE){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->chi[g];
}

inline double mtrl_get_dcoef(const MTRL *m, size_t g)
{
	#ifdef DEBUG
	if(g >= EGSIZE){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->dcoef[g];
}

inline double mtrl_get_sa(const MTRL *m, size_t g)
{
	
	#ifdef DEBUG
	if(g >= EGSIZE){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->sa[g];
}

inline double mtrl_get_sr(const MTRL *m, size_t g)
{
	#ifdef DEBUG
	if(g >= EGSIZE){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->sr[g];
}

inline double mtrl_get_vsf(const MTRL *m, size_t g)
{
	#ifdef DEBUG
	if(g >= EGSIZE){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->vsf[g];
}

inline double mtrl_get_ss(const MTRL *m, size_t g, size_t from_g)
{
	#ifdef DEBUG
	if(g >= EGSIZE){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return m->ss[g][from_g];
}
