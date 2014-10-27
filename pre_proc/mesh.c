#include"mesh.h"
#include<stdlib.h>

double pos_cal(size_t rz_s, size_t rz, const double *rz_span, const size_t *rz_subdiv);
void bdy_check(int ***bdy_checker, int ***mtrl_set, const SCONF *sconf,
		size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t *xspan_subdiv, size_t *yspan_subdiv, size_t *zspan_subdiv);

MESH *mesh_create(const SCONF *sconf)
{
	size_t eg_size = sconf->eg_size;
	size_t xm_mesh_size = sconf->xm_mesh_size;
	size_t ym_mesh_size = sconf->ym_mesh_size;
	size_t zm_mesh_size = sconf->zm_mesh_size;
	MESH *mesh = malloc(sizeof(MESH));
	mesh->eg_size = eg_size;
	mesh->xm_mesh_size = xm_mesh_size;
	mesh->ym_mesh_size = ym_mesh_size;
	mesh->zm_mesh_size = zm_mesh_size;
	// Index by [z,y,x,g,g']
	mesh->mtrl_id = malloc(zm_mesh_size * sizeof(int **));
	mesh->bdy_checker = malloc(zm_mesh_size * sizeof(int **));
	for(size_t k=0; k < zm_mesh_size; ++k){
		mesh->mtrl_id[k] = malloc(ym_mesh_size * sizeof(int *));
		mesh->bdy_checker[k] = malloc(ym_mesh_size * sizeof(int *));
		for(size_t j=0; j < ym_mesh_size; ++j){
			mesh->mtrl_id[k][j] = calloc(xm_mesh_size, sizeof(int));
			mesh->bdy_checker[k][j] = calloc(xm_mesh_size, sizeof(int));
		}
	}
	mesh->dx = cdat3_create(xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->dy = cdat3_create(xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->dz = cdat3_create(xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->xpos = cdat3_create(xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->ypos = cdat3_create(xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->zpos = cdat3_create(xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->chi = cdat4_create(eg_size, xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->dcoef = cdat4_create(eg_size, xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->sa = cdat4_create(eg_size, xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->sr = cdat4_create(eg_size, xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->vsf = cdat4_create(eg_size, xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->ss = cdat5_create(eg_size, xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->adfxl = cdat4_create(eg_size, xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->adfxr = cdat4_create(eg_size, xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->adfyl = cdat4_create(eg_size, xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->adfyr = cdat4_create(eg_size, xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->adfzl = cdat4_create(eg_size, xm_mesh_size, ym_mesh_size, zm_mesh_size);
	mesh->adfzr = cdat4_create(eg_size, xm_mesh_size, ym_mesh_size, zm_mesh_size);
	
	size_t xm_span_size = sconf->xm_span_size;
	size_t ym_span_size = sconf->ym_span_size;
	size_t zm_span_size = sconf->zm_span_size;
	size_t *xspan_subdiv = sconf->xspan_subdiv;
	size_t *yspan_subdiv = sconf->yspan_subdiv;
	size_t *zspan_subdiv = sconf->zspan_subdiv;
	double *xspan_len = sconf->xspan_len;
	double *yspan_len = sconf->yspan_len;
	double *zspan_len = sconf->zspan_len;
	int ***mtrl_set = sconf->mtrl_set;
	MTRLLIB *mlib = sconf->mtrllib;
	bdy_check(mesh->bdy_checker, mtrl_set, sconf, xm_span_size, ym_span_size, zm_span_size, xspan_subdiv, yspan_subdiv, zspan_subdiv);
	for(size_t zspan=0; zspan<zm_span_size; ++zspan)
		for(size_t yspan=0; yspan<ym_span_size; ++yspan)
			for(size_t xspan=0; xspan<xm_span_size; ++xspan){
				MBLOCK mblock = sconf_get_mblock(sconf, xspan, yspan, zspan);
				for(size_t k = mblock.start_z; k <= mblock.end_z; ++k)
					for(size_t j = mblock.start_y; j <= mblock.end_y; ++j)
						for(size_t i = mblock.start_x; i <= mblock.end_x; ++i){
							if(mesh->bdy_checker[k][j][i] & 0b00000001) continue;
							int mtrl_id = mtrl_set[xspan][yspan][zspan];
							MTRL *m = mtrllib_get_fromid(mlib, mtrl_id);
							mesh->mtrl_id[k][j][i] = mtrl_id;
							cdat3_set_val(mesh->dx, i, j, k, xspan_len[xspan]/xspan_subdiv[xspan]);
							cdat3_set_val(mesh->dy, i, j, k, yspan_len[yspan]/yspan_subdiv[yspan]);
							cdat3_set_val(mesh->dz, i, j, k, zspan_len[zspan]/zspan_subdiv[zspan]);
							cdat3_set_val(mesh->xpos, i, j, k, pos_cal(xspan, i, xspan_len, xspan_subdiv));
							cdat3_set_val(mesh->ypos, i, j, k, pos_cal(yspan, j, yspan_len, yspan_subdiv));
							cdat3_set_val(mesh->zpos, i, j, k, pos_cal(zspan, k, zspan_len, zspan_subdiv));
							for(size_t g=0; g< eg_size; ++g){
								cdat4_set_val(mesh->chi, g, i, j, k, m->chi[g]);
								cdat4_set_val(mesh->dcoef, g, i, j, k, m->dcoef[g]);
								cdat4_set_val(mesh->sa, g, i, j, k, m->sa[g]);
								cdat4_set_val(mesh->sr, g, i, j, k, m->sr[g]);
								cdat4_set_val(mesh->vsf, g, i, j, k, m->vsf[g]);
								for(size_t from_g=0; from_g< eg_size; ++from_g)
									cdat5_set_val(mesh->ss, g, from_g, i, j, k, m->ss[g][from_g]);
								// ADFs
								cdat4_set_val(mesh->adfxl, g, i, j, k, (i == mblock.start_x) ? m->adfxl[g] : 1.0);
								cdat4_set_val(mesh->adfxr, g, i, j, k, (i == mblock.end_x) ? m->adfxr[g] : 1.0);
								cdat4_set_val(mesh->adfyl, g, i, j, k, (j == mblock.start_y) ? m->adfyl[g] : 1.0);
								cdat4_set_val(mesh->adfyr, g, i, j, k, (j == mblock.end_y) ? m->adfyr[g] : 1.0);
								cdat4_set_val(mesh->adfzl, g, i, j, k, (k == mblock.start_z) ? m->adfzl[g] : 1.0);
								cdat4_set_val(mesh->adfzr, g, i, j, k, (k == mblock.end_y) ? m->adfzr[g] : 1.0);
							}
						}
			}
	return mesh;
}

void mesh_free(MESH *mesh)
{
	size_t zm_mesh_size = mesh->zm_mesh_size;
	size_t ym_mesh_size = mesh->ym_mesh_size;
	for(size_t k=0; k< zm_mesh_size; ++k){
		for(size_t j=0; j< ym_mesh_size; ++j){
			free(mesh->mtrl_id[k][j]);
			free(mesh->bdy_checker[k][j]);
		}
		free(mesh->mtrl_id[k]);
		free(mesh->bdy_checker[k]);
	}
	free(mesh->mtrl_id);
	free(mesh->bdy_checker);
}

inline int mesh_get_mtrl_id(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->mtrl_id[k][j][i]);
}

inline double mesh_get_dx(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat3_get_val(mesh->dx, i, j, k);
}

inline double mesh_get_dy(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat3_get_val(mesh->dy, i, j, k);
}

inline double mesh_get_dz(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat3_get_val(mesh->dy, i, j, k);
}

inline double mesh_get_xpos(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat3_get_val(mesh->xpos, i, j, k);
}

inline double mesh_get_ypos(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat3_get_val(mesh->ypos, i, j, k);
}

inline double mesh_get_zpos(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat3_get_val(mesh->zpos, i, j, k);
}

inline double mesh_get_chi(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->chi, g, i, j, k);
}

inline double mesh_get_dcoef(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->dcoef, g, i, j, k);
}

inline double mesh_get_sa(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->sa, g, i, j, k);
}

inline double mesh_get_sr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->sr, g, i, j, k);
}

inline double mesh_get_vsf(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->vsf, g, i, j, k);
}

inline double mesh_get_ss(const MESH *mesh, size_t g, size_t from_g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat5_get_val(mesh->ss, g, from_g, i, j, k);
}

inline double mesh_get_adfxl(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->adfxl, g, i, j, k);
}

inline double mesh_get_adfxr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->adfxr, g, i, j, k);
}

inline double mesh_get_adfyl(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->adfyl, g, i, j, k);
}

inline double mesh_get_adfyr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->adfyr, g, i, j, k);
}

inline double mesh_get_adfzl(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->adfzl, g, i, j, k);
}

inline double mesh_get_adfzr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(mesh->bdy_checker[k][j][i] & 0b00000001){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->adfzr, g, i, j, k);
}

void mesh_fprintf(const MESH *mesh, MAPPER *mapper, FILE *stream)
{
	size_t eg_size = mesh->eg_size;
	size_t rt_size = mapper->rt_size;
	fprintf(stream, "1DID\t3DXID\t3DYID\t3DZID\tMTRL_ID\tDX\tDY\tDZ\tXPOS\tYPOS\tZPOS\tXL\tXR\tYL\tYR\tZL\tZR\n");
	for(size_t i=0; i< rt_size; ++i){
		XYZ_IDX xyz = mapper->one2three[i];
		size_t x = xyz.xi;
		size_t y = xyz.yi;
		size_t z = xyz.zi;
		size_t r = mapper->three2one[z][y][x];
		fprintf(stream, "%4zd\t%4zd\t%4zd\t%4zd\t%4d\t%4g\t%4g\t%4g\t%4g\t%4g\t%4g\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\n",
				r, x, y, z, 
				mesh_get_mtrl_id(mesh, x, y, z),
				mesh_get_dx(mesh, x, y, z),
				mesh_get_dy(mesh, x, y, z),
				mesh_get_dz(mesh, x, y, z),
				mesh_get_xpos(mesh, x, y, z),
				mesh_get_ypos(mesh, x, y, z),
				mesh_get_zpos(mesh, x, y, z),
				(mesh->bdy_checker[z][y][x] & 0b00000100) != 0,
				(mesh->bdy_checker[z][y][x] & 0b00001000) != 0,
				(mesh->bdy_checker[z][y][x] & 0b00010000) != 0,
				(mesh->bdy_checker[z][y][x] & 0b00100000) != 0,
				(mesh->bdy_checker[z][y][x] & 0b01000000) != 0,
				(mesh->bdy_checker[z][y][x] & 0b10000000) != 0);
	}
	fprintf(stream, "\n");
	fprintf(stream, "1DID\t4DEID\t4DXID\t4DYID\t4DZID\tADFXL\tADFXR\tADFYL\tADFYR\tADFZL\tADFZR\n");
	for(size_t g=0; g< eg_size; ++g){
		for(size_t i=0; i< rt_size; ++i){
			XYZ_IDX xyz = mapper->one2three[i];
			size_t x = xyz.xi;
			size_t y = xyz.yi;
			size_t z = xyz.zi;
			size_t r = mapper->three2one[z][y][x];
			fprintf(stream, "%4zd\t%4zd\t%4zd\t%4zd\t%4zd\t%4g\t%4g\t%4g\t%4g\t%4g\t%4g\n",
				g*rt_size+r, g, x, y, z,
				mesh_get_adfxl(mesh, g, x, y, z),
				mesh_get_adfxr(mesh, g, x, y, z),
				mesh_get_adfyl(mesh, g, x, y, z),
				mesh_get_adfyr(mesh, g, x, y, z),
				mesh_get_adfzl(mesh, g, x, y, z),
				mesh_get_adfzr(mesh, g, x, y, z));
		}
	}
	fprintf(stream, "\n");
	fprintf(stream, "1DID\t4DEID\t4DXID\t4DYID\t4DZID\tCHI\tDCOEF\tSA\tSR\tVSF\tSS(FROM EG1 ... EGN)\n");
	for(size_t g=0; g< eg_size; ++g){
		for(size_t i=0; i< rt_size; ++i){
			XYZ_IDX xyz = mapper->one2three[i];
			size_t x = xyz.xi;
			size_t y = xyz.yi;
			size_t z = xyz.zi;
			size_t r = mapper->three2one[z][y][x];
			fprintf(stream, "%4zd\t%4zd\t%4zd\t%4zd\t%4zd\t%4.3g\t%4.3g\t%4.3g\t%4.3g\t%4.3g\t",
					g*rt_size+r, g, x, y, z, 
					mesh_get_chi(mesh, g, x, y, z),
					mesh_get_dcoef(mesh, g, x, y, z),
					mesh_get_sa(mesh, g, x, y, z),
					mesh_get_sr(mesh, g, x, y, z),
					mesh_get_vsf(mesh, g, x, y, z));
			for(size_t from_g=0; from_g < eg_size; ++from_g)
				fprintf(stream, "%g\t", mesh_get_ss(mesh, g, from_g, x, y, z));
			fprintf(stream, "\n");
		}
	}
}

// private func.
double pos_cal(size_t rz_s, size_t rz, const double *rz_span, const size_t *rz_subdiv)
{
	double tmp = 0.0;
	double t = rz;
	for(size_t i=0; i< rz_s; ++i){
		tmp += rz_span[i];
		t -= rz_subdiv[i];
	}
	double d = rz_span[rz_s] / rz_subdiv[rz_s];
	tmp += t * d;
	tmp += d/2.0;
	return tmp;
}

// private func.
int has_xlspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan)
{
	if(xspan == 0) return 0;
	if(mtrl_set[xspan-1][yspan][zspan] < 0) return 0;
	return 1;
}

// private func.
int has_xrspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan)
{
	if(xspan == xspan_size-1) return 0;
	if(mtrl_set[xspan+1][yspan][zspan] < 0) return 0;
	return 1;
}

// private func.
int has_ylspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan)
{
	if(yspan == 0) return 0;
	if(mtrl_set[xspan][yspan-1][zspan] < 0) return 0;
	return 1;
}

// private func.
int has_yrspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan)
{
	if(yspan == yspan_size-1) return 0;
	if(mtrl_set[xspan][yspan+1][zspan] < 0) return 0;
	return 1;
}

// private func.
int has_zlspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan)
{
	if(zspan == 0) return 0;
	if(mtrl_set[xspan][yspan][zspan-1] < 0) return 0;
	return 1;
}

// private func.
int has_zrspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan)
{
	if(zspan == zspan_size-1) return 0;
	if(mtrl_set[xspan][yspan][zspan+1] < 0) return 0;
	return 1;
}

// private func.
void bdy_check(int ***bdy_checker, int ***mtrl_set, const SCONF *sconf,
		size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t *xspan_subdiv, size_t *yspan_subdiv, size_t *zspan_subdiv)
{
	for(size_t zspan=0; zspan<zspan_size; ++zspan)
		for(size_t yspan=0; yspan<yspan_size; ++yspan)
			for(size_t xspan=0; xspan<xspan_size; ++xspan){
				MBLOCK mblock = sconf_get_mblock(sconf, xspan, yspan, zspan);
				for(size_t k = mblock.start_z; k <= mblock.end_z; ++k)
					for(size_t j = mblock.start_y; j <= mblock.end_y; ++j)
						for(size_t i = mblock.start_x; i <= mblock.end_x; ++i){
						int mtrl_id = mtrl_set[xspan][yspan][zspan];
						if(mtrl_id < 0){
							bdy_checker[k][j][i] = 0b00000001;
							continue;
						}
						else
							bdy_checker[k][j][i] = 0b00000010;

						if(i == mblock.start_x && !has_xlspan(mtrl_set,xspan_size,yspan_size,zspan_size,xspan,yspan,zspan))
							bdy_checker[k][j][i] |= 0b00000100;
						if(i == mblock.end_x && !has_xrspan(mtrl_set,xspan_size,yspan_size,zspan_size,xspan,yspan,zspan))
							bdy_checker[k][j][i] |= 0b00001000;
						if(j == mblock.start_y && !has_ylspan(mtrl_set,xspan_size,yspan_size,zspan_size,xspan,yspan,zspan))
							bdy_checker[k][j][i] |= 0b00010000;
						if(j == mblock.end_y && !has_yrspan(mtrl_set,xspan_size,yspan_size,zspan_size,xspan,yspan,zspan))
							bdy_checker[k][j][i] |= 0b00100000;
						if(k == mblock.start_z && !has_zlspan(mtrl_set,xspan_size,yspan_size,zspan_size,xspan,yspan,zspan))
							bdy_checker[k][j][i] |= 0b01000000;
						if(k == mblock.end_z && !has_zrspan(mtrl_set,xspan_size,yspan_size,zspan_size,xspan,yspan,zspan))
							bdy_checker[k][j][i] |= 0b10000000;
						}
			}
}
