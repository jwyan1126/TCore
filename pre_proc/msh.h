#ifndef MSH_H
#define MSH_H

#include"mapper.h"
#include<stddef.h>
#include"sconf.h"

typedef struct
{
	int *mtrl_id;
	double *chi;
	double *dcoef;
	double *sa;
	double *sr;
	double *vsf;
	double **ss;
	double *dx;
	double *dy;
	double *dz;
	double *xpos;
	double *ypos;
	double *zpos;
	size_t eg_size;
	size_t rt_size;
	MAPPER *mapper;
} MSH;

MSH *msh_create(const SCONF *sconf, MAPPER *mapper);

void msh_free(MSH *msh);

int msh_get_mtrl_id(const MSH *msh, size_t i, size_t j, size_t k);

double msh_get_chi(const MSH *msh, size_t g, size_t i, size_t j, size_t k);

double msh_get_dcoef(const MSH *msh, size_t g, size_t i, size_t j, size_t k);

double msh_get_sa(const MSH *msh, size_t g, size_t i, size_t j, size_t k);

double msh_get_sr(const MSH *msh, size_t g, size_t i, size_t j, size_t k);

double msh_get_vsf(const MSH *msh, size_t g, size_t i, size_t j, size_t k);

double msh_get_ss(const MSH *msh, size_t g, size_t from_g, size_t i, size_t j, size_t k);

double msh_get_dx(const MSH *msh, size_t i, size_t j, size_t k);

double msh_get_dy(const MSH *msh, size_t i, size_t j, size_t k);

double msh_get_dz(const MSH *msh, size_t i, size_t j, size_t k);

double msh_get_xpos(const MSH *msh, size_t i, size_t j, size_t k);

double msh_get_ypos(const MSH *msh, size_t i, size_t j, size_t k);

double msh_get_zpos(const MSH *msh, size_t i, size_t j, size_t k);

void msh_fprintf(const MSH *msh, FILE *stream);

#endif
