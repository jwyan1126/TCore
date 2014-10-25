#include"mesh.h"

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
	mesh->dx = malloc(zm_mesh_size * sizeof(double **));
	mesh->dy = malloc(zm_mesh_size * sizeof(double **));
	mesh->dz = malloc(zm_mesh_size * sizeof(double **));
	mesh->xpos = malloc(zm_mesh_size * sizeof(double **));
	mesh->ypos = malloc(zm_mesh_size * sizeof(double **));
	mesh->zpos = malloc(zm_mesh_size * sizeof(double **));
	mesh->chi = malloc(zm_mesh_size * sizeof(double ***));
	mesh->dcoef = malloc(zm_mesh_size * sizeof(double ***));
	mesh->sa = malloc(zm_mesh_size * sizeof(double ***));
	mesh->sr = malloc(zm_mesh_size * sizeof(double ***));
	mesh->vsf = malloc(zm_mesh_size * sizeof(double ***));
	mesh->ss = malloc(zm_mesh_size * sizeof(double ****));
	mesh->adfxl = malloc(zm_mesh_size * sizeof(double ***));
	mesh->adfxr = malloc(zm_mesh_size * sizeof(double ***));
	mesh->adfyl = malloc(zm_mesh_size * sizeof(double ***));
	mesh->adfyr = malloc(zm_mesh_size * sizeof(double ***));
	mesh->adfzl = malloc(zm_mesh_size * sizeof(double ***));
	mesh->adfzr = malloc(zm_mesh_size * sizeof(double ***));
	for(size_t k=0; k < zm_mesh_size; ++k){
		mesh->mtrl_id[k] = malloc(ym_mesh_size * sizeof(int *));
		mesh->check[k] = malloc(ym_mesh_size * sizeof(int *));
		mesh->dx[k] = malloc(ym_mesh_size * sizeof(double *));
		mesh->dy[k] = malloc(ym_mesh_size * sizeof(double *));
		mesh->dz[k] = malloc(ym_mesh_size * sizeof(double *));
		mesh->xpos[k] = malloc(ym_mesh_size * sizeof(double *));
		mesh->ypos[k] = malloc(ym_mesh_size * sizeof(double *));
		mesh->zpos[k] = malloc(ym_mesh_size * sizeof(double *));
		mesh->chi[k] = malloc(ym_mesh_size * sizeof(double **));
		mesh->dcoef[k] = malloc(ym_mesh_size * sizeof(double **));
		mesh->sa[k] = malloc(ym_mesh_size * sizeof(double **));
		mesh->sr[k] = malloc(ym_mesh_size * sizeof(double **));
		mesh->vsf[k] = malloc(ym_mesh_size * sizeof(double **));
		mesh->ss[k] = malloc(ym_mesh_size * sizeof(double ***));
		mesh->adfxl[k] = malloc(ym_mesh_size * sizeof(double **));
		mesh->adfxr[k] = malloc(ym_mesh_size * sizeof(double **));
		mesh->adfyl[k] = malloc(ym_mesh_size * sizeof(double **));
		mesh->adfyr[k] = malloc(ym_mesh_size * sizeof(double **));
		mesh->adfzl[k] = malloc(ym_mesh_size * sizeof(double **));
		mesh->adfzr[k] = malloc(ym_mesh_size * sizeof(double **));
		for(size_t j=0; j < ym_mesh_size; ++j){
			mesh->mtrl_id[k][j] = calloc(xm_mesh_size, sizeof(int));
			mesh->check[k][j] = calloc(xm_mesh_size, sizeof(int));
			mesh->dx[k][j] = calloc(xm_mesh_size, sizeof(double));
			mesh->dy[k][j] = calloc(xm_mesh_size, sizeof(double));
			mesh->dz[k][j] = calloc(xm_mesh_size, sizeof(double));
			mesh->xpos[k][j] = calloc(xm_mesh_size, sizeof(double));
			mesh->ypos[k][j] = calloc(xm_mesh_size, sizeof(double));
			mesh->zpos[k][j] = calloc(xm_mesh_size, sizeof(double));
			mesh->chi[k][j] = malloc(xm_mesh_size * sizeof(double *));
			mesh->dcoef[k][j] = malloc(xm_mesh_size * sizeof(double *));
			mesh->sa[k][j] = malloc(xm_mesh_size * sizeof(double *));
			mesh->sr[k][j] = malloc(xm_mesh_size * sizeof(double *));
			mesh->vsf[k][j] = malloc(xm_mesh_size * sizeof(double *));
			mesh->ss[k][j] = malloc(xm_mesh_size * sizeof(double **));
			mesh->adfxl[k][j] = malloc(xm_mesh_size * sizeof(double *));
			mesh->adfxr[k][j] = malloc(xm_mesh_size * sizeof(double *));
			mesh->adfyl[k][j] = malloc(xm_mesh_size * sizeof(double *));
			mesh->adfyr[k][j] = malloc(xm_mesh_size * sizeof(double *));
			mesh->adfzl[k][j] = malloc(xm_mesh_size * sizeof(double *));
			mesh->adfzr[k][j] = malloc(xm_mesh_size * sizeof(double *));
			for(size_t g=0; g< eg_size; ++g){
				mesh->chi[k][j][i] = calloc(eg_size, sizeof(double));
				mesh->dcoef[k][j][i] = calloc(eg_size, sizeof(double));
				mesh->sa[k][j][i] = calloc(eg_size, sizeof(double));
				mesh->sr[k][j][i] = calloc(eg_size, sizeof(double));
				mesh->vsf[k][j][i] = calloc(eg_size, sizeof(double));
				mesh->ss[k][j][i] = malloc(eg_size * sizeof(double *));
				mesh->adfxl[k][j][i] = calloc(eg_size, sizeof(double));
				mesh->adfxr[k][j][i] = calloc(eg_size, sizeof(double));
				mesh->adfyl[k][j][i] = calloc(eg_size, sizeof(double));
				mesh->adfyr[k][j][i] = calloc(eg_size, sizeof(double));
				mesh->adfzl[k][j][i] = calloc(eg_size, sizeof(double));
				mesh->adfzr[k][j][i] = calloc(eg_size, sizeof(double));
				for(size_t from_g=0; from_g< eg_size; ++from_g)
					mesh->ss[k][j][i][g] = calloc(eg_size, sizeof(double));
			}
		}
	}
	
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
	for(size_t k=0; k<zm_span_size; ++k)
		for(size_t j=0; j<ym_span_size; ++j)
			for(size_t i=0; i<xm_span_size; ++i){
				MBLOCK //////
			}
	return mesh;
}

void mesh_free(MESH *mesh)
{
	for(size_t k=0; k< zm_mesh_size; ++k){
		for(size_t j=0; j< ym_mesh_size; ++j){
			for(size_t i=0; i< xm_mesh_size; ++i){
				for(size_t g=0; g< eg_size; ++g)
					free(mesh->ss[k][j][i][g]);
				free(mesh->chi[k][j][i]);
				free(mesh->dcoef[k][j][i]);
				free(mesh->sa[k][j][i]);
				free(mesh->sr[k][j][i]);
				free(mesh->vsf[k][j][i]);
				free(mesh->ss[k][j][i]);
				free(mesh->adfxl[k][j][i]);
				free(mesh->adfxr[k][j][i]);
				free(mesh->adfyl[k][j][i]);
				free(mesh->adfyr[k][j][i]);
				free(mesh->adfzl[k][j][i]);
				free(mesh->adfzr[k][j][i]);
			}
			free(mesh->mtrl_id[k][j]);
			free(mesh->check[k][j]);
			free(mesh->dx[k][j]);
			free(mesh->dy[k][j]);
			free(mesh->dz[k][j]);
			free(mesh->xpos[k][j]);
			free(mesh->ypos[k][j]);
			free(mesh->zpos[k][j]);
			free(mesh->chi[k][j]);
			free(mesh->dcoef[k][j]);
			free(mesh->sa[k][j]);
			free(mesh->sr[k][j]);
			free(mesh->vsf[k][j]);
			free(mesh->ss[k][j]);
			free(mesh->adfxl[k][j]);
			free(mesh->adfxr[k][j]);
			free(mesh->adfyl[k][j]);
			free(mesh->adfyr[k][j]);
			free(mesh->adfzl[k][j]);
			free(mesh->adfzr[k][j]);
		}
		free(mesh->mtrl_id[k]);
		free(mesh->dx[k]);
		free(mesh->dy[k]);
		free(mesh->dz[k]);
		free(mesh->xpos[k]);
		free(mesh->ypos[k]);
		free(mesh->zpos[k]);
		free(mesh->chi[k]);
		free(mesh->dcoef[k]);
		free(mesh->sa[k]);
		free(mesh->sr[k]);
		free(mesh->vsf[k]);
		free(mesh->ss[k]);
		free(mesh->adfxl[k]);
		free(mesh->adfxr[k]);
		free(mesh->adfyl[k]);
		free(mesh->adfyr[k]);
		free(mesh->adfzl[k]);
		free(mesh->adfzr[k]);
	}
	free(mesh->mtrl_id);
	free(mesh->dx);
	free(mesh->dy);
	free(mesh->dz);
	free(mesh->xpos);
	free(mesh->ypos);
	free(mesh->zpos);
	free(mesh->chi);
	free(mesh->dcoef);
	free(mesh->sa);
	free(mesh->sr);
	free(mesh->vsf);
	free(mesh->ss);
	free(mesh->adfxl);
	free(mesh->adfxr);
	free(mesh->adfyl);
	free(mesh->adfyr);
	free(mesh->adfzl);
	free(mesh->adfzr);
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
	return(mesh->dx[k][j][i]);
}

inline double mesh_get_dy(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->dy[k][j][i]);
}

inline double mesh_get_dz(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->dz[k][j][i]);
}

inline double mesh_get_xpos(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->xpos[k][j][i]);
}

inline double mesh_get_ypos(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->ypos[k][j][i]);
}

inline double mesh_get_zpos(const MESH *mesh, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->zpos[k][j][i]);
}

inline double mesh_get_chi(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->chi[k][j][i][g]);
}

inline double mesh_get_dcoef(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->dcoef[k][j][i][g]);
}

inline double mesh_get_sa(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->sa[k][j][i][g]);
}

inline double mesh_get_sr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->sr[k][j][i][g]);
}

inline double mesh_get_vsf(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->vsf[k][j][i][g]);
}

inline double mesh_get_ss(const MESH *mesh, size_t g, size_t from_g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->ss[k][j][i][g][from_g]);
}

inline double mesh_get_adfxl(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->adfxl[k][j][i][g]);
}

inline double mesh_get_adfxr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->adfxr[k][j][i][g]);
}

inline double mesh_get_adfyl(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->adfyl[k][j][i][g]);
}

inline double mesh_get_adfyr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->adfyr[k][j][i][g]);
}

inline double mesh_get_adfzl(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->adfzl[k][j][i][g]);
}

inline double mesh_get_adfzr(const MESH *mesh, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(!mesh->check[k][j][i]){
		fprintf(stderr, "Specified mesh is not filled in.\n");
		exit(-1);
	}
	#endif
	return(mesh->adfzr[k][j][i][g]);
}
