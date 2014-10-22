#ifndef CONF_H
#define CONF_H

typedef struct
{
	size_t eg_size; // num of energy groups
	size_t xm_span_size; // max num of spans in x direction
	size_t ym_span_size; // max num of spans in y direction
	size_t zm_span_size; // max num of spans in z direction
	size_t xm_mesh_size; // max num of meshes in x direction
	size_t ym_mesh_size; // max num of meshes in y direction
	size_t zm_mesh_size; // max num of meshes in z direction
	double *xspan_len; // each length of span in x direction
	double *yspan_len; // each length of span in y direction
	double *zspan_len; // each length of span in z direction
	size_t *xspan_subdiv; // num of sub-division of each span in x direction
	size_t *yspan_subdiv; // num of sub-division of each span in y direction
	size_t *zspan_subdiv; // num of sub-division of each span in z direction

	double *xl_bdy; // albedos of each energy group in XL boundary
	double *xr_bdy; // albedos of each energy group in XR boundary
	double *yl_bdy; // albedos of each energy group in YL boundary
	double *yr_bdy; // albedos of each energy group in YR boundary
	double *zl_bdy; // albedos of each energy group in ZL boundary
	double *zr_bdy; // albedos of each energy group in ZR boundary
	
	int ***mtrl_set; // setting of materials in each spanned block
	MTRLLIB *mtrllib; // material library
	
} SCONF; //steady configuration struct

SCONF *sconf_create(INPUT *input);

void sconf_free(SCONF *sconf);

#endif
