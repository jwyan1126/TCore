#ifndef RECT_MSH_H
#define RECT_MSH_H

typedef struct
{
	double *chi_mesh;
	double *dcoef_mesh;
	double *sa_mesh;
	double *vsf_mesh;
	double **ss_mesh;
	double *dx_mesh;
	double *dy_mesh;
	double *dz_mesh;
	double *xpos_mesh;
	double *ypos_mesh;
	double *zpos_mesh;
} RECT_MSH; 

RECT_MSH *rectmsh_create(RECT_MAPPER *rect_mapper/*OUT*/, SCONF *sconf/*IN*/);

void rectmsh_free(RECTMSH *rectmsh);

#endif
