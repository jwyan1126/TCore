#ifndef MSH_H
#define MSH_H

#include"rect_mapper.h"
#include<stddef.h>
#include"sconf.h"

//extern size_t EG_SIZE;
//extern size_t XM_SIZE;
//extern size_t YM_SIZE;
//extern size_t ZM_SIZE;
//extern size_t RT_SIZE;
//extern RECT_MAPPER *MAPPER;

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
} MSH;

MSH *msh_create(const SCONF *sconf);

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

#endif
