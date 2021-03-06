#ifndef MESH_H
#define MESH_H

#include<stddef.h>
#include"sconf.h"
#include"cdat.h"
#include"mapper.h"

typedef struct
{
	size_t eg_size;
	size_t xm_mesh_size;
	size_t ym_mesh_size;
	size_t zm_mesh_size;
	size_t rt_size;
	MAPPER *mapper;
	int ***cchecker;
	int ***xchecker;
	int ***ychecker;
	int ***zchecker;
	int ***mtrl_id; // indexed by [z][y][x]
	int xl_bdy;
	int xr_bdy;
	int yl_bdy;
	int yr_bdy;
	int zl_bdy;
	int zr_bdy;
	CDAT3 *dx;
	CDAT3 *dy;
	CDAT3 *dz;
	CDAT3 *xpos;
	CDAT3 *ypos;
	CDAT3 *zpos;
	CDAT4 *chi;
	CDAT4 *dcoef;
	CDAT4 *sa;
	CDAT4 *sr;
	CDAT4 *vsf;
	CDAT5 *ss;
	CDAT4 *adfxl;
	CDAT4 *adfxr;
	CDAT4 *adfyl;
	CDAT4 *adfyr;
	CDAT4 *adfzl;
	CDAT4 *adfzr;
} MESH;

MESH *mesh_create(SCONF *sconf, MAPPER *mapper);

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

void mesh_fprintf(const MESH *mesh, FILE *stream);
#endif
