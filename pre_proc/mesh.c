#include"mesh.h"

double pos_cal(size_t rz_s, size_t rz, const double *rz_span, const size_t *rz_subdiv);

MESH *mesh_create(const SCONF *sconf)
{
	size_t eg_size = sconf->eg_size;
	size_t xm_mesh_size = sconf->xm_mesh_size;
	size_t ym_mesh_size = sconf->ym_mesh_size;
	size_t zm_mesh_size = sconf->zm_mesh_size;
	MESH mesh = malloc(sizeof(MESH));
	mesh->eg_size = eg_size;
	mesh->xm_mesh_size = xm_mesh_size;
	mesh->ym_mesh_size = ym_mesh_size;
	mesh->zm_mesh_size = zm_mesh_size;
	// Index by [z,y,x,g,g']
	mesh->mtrl_id = malloc(zm_mesh_size * sizeof(int **));
	mesh->check = malloc(zm_mesh_size * sizeof(int **));
	for(size_t k=0; k < zm_mesh_size; ++k){
		mesh->mtrl_id[k] = malloc(ym_mesh_size * sizeof(int *));
		mesh->check[k] = malloc(ym_mesh_size * sizeof(int *));
		for(size_t j=0; j < ym_mesh_size; ++j){
			mesh->mtrl_id[k][j] = calloc(xm_mesh_size, sizeof(int));
			mesh->check[k][j] = calloc(xm_mesh_size, sizeof(int));
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
	for(size_t zspan=0; zspan<zm_span_size; ++zspan)
		for(size_t yspan=0; yspan<ym_span_size; ++yspan)
			for(size_t xspan=0; xspan<xm_span_size; ++xspan){
				MBLOCK mblock = sconf_get_mblock(sconf, xspan, yspan, zspan);
				for(size_t k = mblock.start_z; k <= mblock.end_z; ++k)
					for(size_t j = mblock.start_y; j <= mblock.end_y; ++j)
						for(size_t i = mblock.start_x; i <= mblock.end_x; ++i){
							mesh->check[k][j][i] = 1;
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
	for(size_t k=0; k< zm_mesh_size; ++k){
		for(size_t j=0; j< ym_mesh_size; ++j){
			free(mesh->mtrl_id[k][j]);
			free(mesh->check[k][j]);
		}
		free(mesh->mtrl_id[k]);
		free(mesh->check[k]);
	}
	free(mesh->mtrl_id);
	free(mesh->check);
}

inline int mesh_get_mtrl_id(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->mtrl_id[k][j][i]);
}

inline double mesh_get_dx(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat3_get_val(mesh->dx, i, j, k);
}

inline double mesh_get_dy(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat3_get_val(mesh->dy, i, j, k);
}

inline double mesh_get_dz(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat3_get_val(mesh->dy, i, j, k);
}

inline double mesh_get_xpos(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat3_get_val(mesh->xpos, i, j, k);
}

inline double mesh_get_ypos(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat3_get_val(mesh->ypos, i, j, k);
}

inline double mesh_get_zpos(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat3_get_val(mesh->zpos, i, j, k);
}

inline double mesh_get_chi(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->chi, g, i, j, k);
}

inline double mesh_get_dcoef(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->dcoef, g, i, j, k);
}

inline double mesh_get_sa(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->sa, g, i, j, k);
}

inline double mesh_get_sr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->sr, g, i, j, k);
}

inline double mesh_get_vsf(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->vsf, g, i, j, k);
}

inline double mesh_get_ss(const MESH *mesh, size_t g, size_t from_g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat5_get_val(mesh->ss, g, from_g, i, j, k);
}

inline double mesh_get_adfxl(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->adfxl, g, i, j, k);
}

inline double mesh_get_adfxr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->adfxr, g, i, j, k);
}

inline double mesh_get_adfyl(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->adfyl, g, i, j, k);
}

inline double mesh_get_adfyr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->adfyr, g, i, j, k);
}

inline double mesh_get_adfzl(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->adfzl, g, i, j, k);
}

inline double mesh_get_adfzr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return cdat4_get_val(mesh->adfzr, g, i, j, k);
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
