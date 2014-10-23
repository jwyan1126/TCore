#include"msh.h"
#include"rect_mapper.h"

size_t inwhichspan(size_t arr[], size_t span_len, int i);
double rzpos_cal(size_t rz_s/*第几个span*/, size_t rz/*第几个网格*/, const double *rz_span/*各span长度*/, const size_t *rz_subdiv/*各span网格数*/);

MSH *msh_create(const SCONF *sconf)
{
	MSH *msh = malloc(sizeof(MSH));
	msh->mtrl_id = malloc(RT_SIZE * sizeof(int));
	msh->chi = malloc(EG_SIZE * RT_SIZE * sizeof(double));
	msh->dcoef = malloc(EG_SIZE * RT_SIZE * sizeof(double));
	msh->sa = malloc(EG_SIZE * RT_SIZE * sizeof(double));
	msh->sr = malloc(EG_SIZE * RT_SIZE * sizeof(double));
	msh->ss = malloc(EG_SIZE * RT_SIZE * sizeof(double *));
	for(size_t i=0; i<EG_SIZE * RT_SIZE; ++i)
		msh->ss[i] = malloc(EG_SIZE & sizeof(double));
	msh->dx = malloc(RT_SIZE * sizeof(double));
	msh->dy = malloc(RT_SIZE * sizeof(double));
	msh->dz = malloc(RT_SIZE * sizeof(double));
	msh->xpos = malloc(RT_SIZE * sizeof(double));
	msh->ypos = malloc(RT_SIZE * sizeof(double));
	msh->zpos = malloc(RT_SIZE * sizeof(double));
	
	// traversal
	size_t xm_span_size = sconf->xm_span_size;
	size_t ym_span_size = sconf->ym_span_size;
	size_t zm_span_size = sconf->zm_span_size;
	size_t xm_mesh_size = sconf->xm_mesh_size;
	size_t ym_mesh_size = sconf->ym_mesh_size;
	size_t zm_mesh_size = sconf->zm_mesh_size;
	size_t ac = 0; // array counter
	MTRLLIB *mlib = sconf->mtrllib;
	for(size_t k = 0; k < zm_mesh_size; ++k){
		for(size_t j = 0; j < ym_mesh_size; ++j){
			for(size_t i = 0; i < xm_mesh_size; ++i){
				size_t xspan = inwhichspan(sconf->xspan_subdiv, sconf->xm_span_size, i);
				size_t yspan = inwhichspan(sconf->yspan_subdiv, sconf->ym_span_size, j);
				size_t zspan = inwhichspan(sconf->zspan_subdiv, sconf->zm_span_size, k);
				int mtrl_id = sconf->mtrl_set[i][j][k];
				if(mtrl_id < 0) continue;
				MTRL *m = mtrllib_get_fromid(mlib, mtrl_id);
				MAPPER->three2one[k][j][i] = ac;
				XYZ_IDX xyz;
				xyz.xi = i; xyz.yi = j; xyz.zi = k;
				MAPPER->one2three[ac] = xyz;
				msh->mtrl_id[ac] = mtrl_id;
				// dx dy dz xpos ypos zpos
				// ......
				//
				for(size_t g=0; g<EG_SIZE; ++g){
					msh->chi[g*RT_SIZE+ac] = m->chi[g];
					msh->dcoef[g*RT_SIZE+ac] = m->dcoef[g];
					msh->sa[g*RT_SIZE+ac] = m->sa[g];
					msh->sr[g*RT_SIZE+ac] = m->sr[g];
					msh->vsf[g*RT_SIZE+ac] = m->vsf[g];
					for(size_t from_g=0; from_g<EG_SIZE; ++from_g)
						msh->ss[g*RT_SIZE+ac][from_g] = m->ss[g][from_g];
				}
				++ac;
			}
		}
	}
	return msh;
}

void msh_free(MSH msh)
{
	free(msh->mtrl_id);
	free(msh->chi);
	free(msh->dcoef);
	free(msh->sa);
	free(msh->vsf);
	for(size_t i=0; i<RT_SIZE; ++i)
		free(msh->ss[i]);
	free(msh->ss);
	free(msh->dx);
	free(msh->dy);
	free(msh->dz);
	free(msh->xpos);
	free(msh->ypos);
	free(msh->zpos);
	free(msh);
}

inline int msh_get_mtrl_id(const MSH *msh, size_t i, size_t j, size_t k)
{
	return msh->mtrl_id[rect_mapper_get1Didx(MAPPER, i,j,k)];
}

inline double msh_get_chi(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->chi[g*RT_SIZE + rect_mapper_get1Didx(MAPPER,i,j,k)];
}

inline double msh_get_dcoef(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->dcoef[g*RT_SIZE + rect_mapper_get1Didx(MAPPER,i,j,k)];
}

inline double msh_get_sa(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->sa[g*RT_SIZE + rect_mapper_get1Didx(MAPPER,i,j,k)];
}

inline double msh_get_sr(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->sr[g*RT_SIZE + rect_mapper_get1Didx(MAPPER,i,j,k)];
}

inline double msh_get_vsf(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->vsf[g*RT_SIZE + rect_mapper_get1Didx(MAPPER,i,j,k)];
}

inline double msh_get_ss(const MSH *msh, size_t g, size_t from_g, size_t i, size_t j, size_t k)
{
	return msh->ss[g*RT_SIZE + rect_mapper_get1Didx(MAPPER,i,j,k)][from_g];
}

inline double msh_get_dx(const MSH *msh, size_t i, size_t j, size_t k)
{
	return msh->dx[g*RT_SIZE + rect_mapper_get1Didx(MAPPER,i,j,k)];
}

inline double msh_get_dy(const MSH *msh, size_t i, size_t j, size_t k)
{
	return msh->dy[g*RT_SIZE + rect_mapper_get1Didx(MAPPER,i,j,k)];
}

inline double msh_get_dz(const MSH *msh, size_t i, size_t j, size_t k)
{
	return msh->dz[g*RT_SIZE + rect_mapper_get1Didx(MAPPER,i,j,k)];
}

inline double msh_get_xpos(const MSH *msh, size_t i, size_t j, size_t k)
{
	return msh->xpos[g*RT_SIZE + rect_mapper_get1Didx(MAPPER,i,j,k)];
}

inline double msh_get_ypos(const MSH *msh, size_t i, size_t j, size_t k)
{
	return msh->ypos[g*RT_SIZE + rect_mapper_get1Didx(MAPPER,i,j,k)];
}

inline double msh_get_zpos(const MSH *msh, size_t i, size_t j, size_t k)
{
	return msh->zpos[g*RT_SIZE + rect_mapper_get1Didx(MAPPER,i,j,k)];
}

// private func.
size_t inwhichspan(size_t arr[], size_t span_len, int i)
{
	for(size_t k=0; k<span_len; ++k){
		i -= arr[k];
		if(i < 0)
			return k;
	}
	fprintf(stderr, "Error occurs. when invoke 'inwhichspan'.\n");
	exit(-1);
}
// private func.
double rzpos_cal(size_t rz_s/*第几个span*/, size_t rz/*第几个网格*/, const double *rz_span/*各span长度*/, const size_t *rz_subdiv/*各span网格数*/)
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
