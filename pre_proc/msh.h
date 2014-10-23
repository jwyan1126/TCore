#ifndef MSH_H
#define MSH_H

extern size_t EG_SIZE;
extern size_t XM_SIZE;
extern size_t YM_SIZE;
extern size_t ZM_SIZE;
extern size_t RT_SIZE;
extern RECT_MAPPER MAPPER;

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

MSH *msh_create(SCONF sconf);

void msh_free(MSH msh);

int msh_get_mtrl_id(size_t i, size_t j, size_t k);

double msh_get_chi(size_t g, size_t i, size_t j, size_t k);

double msh_get_dcoef(size_t g, size_t i, size_t j, size_t k);

double msh_get_sa(size_t g, size_t i, size_t j, size_t k);

double msh_get_sr(size_t g, size_t i, size_t j, size_t k);

double msh_get_vsf(size_t g, size_t i, size_t j, size_t k);

double msh_get_ss(size_t g, size_t from_g, size_t i, size_t j, size_t k);

double msh_get_dx(size_t i, size_t j, size_t k);

double msh_get_dy(size_t i, size_t j, size_t k);

double msh_get_dz(size_t i, size_t j, size_t k);

double msh_get_xpos(size_t i, size_t j, size_t k);

double msh_get_ypos(size_t i, size_t j, size_t k);

double msh_get_zpos(size_t i, size_t j, size_t k);

#endif
