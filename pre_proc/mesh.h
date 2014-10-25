#ifndef MESH_H
#define MESH_H

#include<stddef.h>
#include<sconf.h>

typedef struct
{
	size_t eg_size;
	size_t xm_size;
	size_t ym_size;
	size_t zm_size;
	int ***check;
	int ***mtrl_id;
	double ***dx;
	double ***dy;
	double ***dz;
	double ***xpos;
	double ***ypos;
	double ***zpos;
	double ****chi;
	double ****dcoef;
	double ****sa;
	double ****sr;
	double ****vsf;
	double *****ss;
	double ****adfxl;
	double ****adfxr;
	double ****adfyl;
	double ****adfyr;
	double ****adfzl;
	double ****adfzr;
} MESH;

MESH *mesh_create(const SCONF *sconf);

void mesh_free(MESH *mesh);

int mesh_get_mtrl_id(const MESH *mesh, size_t i, size_t j, size_t k);

double mesh_get_dx(const MESH *mesh, size_t i, size_t j, size_t k);

double mesh_get_dy(const MESH *mesh, size_t i, size_t j, size_t k);

double mesh_get_dz(const MESH *mesh, size_t i, size_t j, size_t k);

double mesh_get_xpos(const MESH *mesh, size_t i, size_t j, size_t k);

double mesh_get_ypos(const MESH *mesh, size_t i, size_t j, size_t k);

double mesh_get_zpos(const MESH *mesh, size_t i, size_t j, size_t k);

double mesh_get_chi(const MESH *mesh, size_t g, size_t i, size_t j, size_t k);

double mesh_get_dcoef(const MESH *mesh, size_t g, size_t i, size_t j, size_t k);

double mesh_get_sa(const MESH *mesh, size_t g, size_t i, size_t j, size_t k);

double mesh_get_sr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k);

double mesh_get_vsf(const MESH *mesh, size_t g, size_t i, size_t j, size_t k);

double mesh_get_ss(const MESH *mesh, size_t g, size_t from_g, size_t i, size_t j, size_t k);

double mesh_get_adfxl(const MESH *mesh, size_t g, size_t i, size_t j, size_t k);

double mesh_get_adfxr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k);

double mesh_get_adfyl(const MESH *mesh, size_t g, size_t i, size_t j, size_t k);

double mesh_get_adfyr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k);

double mesh_get_adfzl(const MESH *mesh, size_t g, size_t i, size_t j, size_t k);

double mesh_get_adfzr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k);
#endif
