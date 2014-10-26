#include"msh.h"
#include"mapper.h"
#include<stdio.h>
#include<stdlib.h>
#include"mtrllib.h"
#include"mtrl.h"

size_t inwhichspan(size_t arr[], size_t span_len, int i);
double pos_cal(size_t rz_s, size_t rz, const double *rz_span, const size_t *rz_subdiv);

MSH *msh_create(const SCONF *sconf, MAPPER *mapper)
{
	size_t eg_size = sconf->eg_size;
	size_t rt_size = sconf->rt_mesh_size;
	MSH *msh = malloc(sizeof(MSH));
	msh->rt_size = rt_size;
	msh->eg_size = eg_size;
	msh->mapper = mapper;
	msh->mtrl_id = malloc(rt_size * sizeof(int));
	msh->dx = malloc(rt_size * sizeof(double));
	msh->dy = malloc(rt_size * sizeof(double));
	msh->dz = malloc(rt_size * sizeof(double));
	msh->xpos = malloc(rt_size * sizeof(double));
	msh->ypos = malloc(rt_size * sizeof(double));
	msh->zpos = malloc(rt_size * sizeof(double));
	msh->chi = malloc(eg_size * rt_size * sizeof(double));
	msh->dcoef = malloc(eg_size * rt_size * sizeof(double));
	msh->sa = malloc(eg_size * rt_size * sizeof(double));
	msh->sr = malloc(eg_size * rt_size * sizeof(double));
	msh->vsf = malloc(eg_size * rt_size * sizeof(double));
	msh->ss = malloc(eg_size * rt_size * sizeof(double *));
	msh->adfxl = malloc(eg_size * rt_size * sizeof(double));
	msh->adfxr = malloc(eg_size * rt_size * sizeof(double));
	msh->adfyl = malloc(eg_size * rt_size * sizeof(double));
	msh->adfyr = malloc(eg_size * rt_size * sizeof(double));
	msh->adfzl = malloc(eg_size * rt_size * sizeof(double));
	msh->adfzr = malloc(eg_size * rt_size * sizeof(double));
	for(size_t i=0; i<eg_size * rt_size; ++i)
		msh->ss[i] = malloc(eg_size * sizeof(double));
	
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
	for(size_t i=0; i < rt_size; ++i){
		XYZ_IDX xyz = mapper_get3Didx(mapper, i);
		size_t x = xyz.xi;
		size_t y = xyz.yi;
		size_t z = xyz.zi;
		size_t xspan = inwhichspan(xspan_subdiv, xm_span_size, x);
		size_t yspan = inwhichspan(yspan_subdiv, ym_span_size, y);
		size_t zspan = inwhichspan(zspan_subdiv, zm_span_size, z);
		int mtrl_id = mtrl_set[xspan][yspan][zspan];
		MTRL *m = mtrllib_get_fromid(mlib, mtrl_id);
		msh->mtrl_id[i] = mtrl_id;
		msh->dx[i] = xspan_len[xspan] / xspan_subdiv[xspan];
		msh->dy[i] = yspan_len[yspan] / yspan_subdiv[yspan];
		msh->dz[i] = zspan_len[zspan] / zspan_subdiv[zspan];
		msh->xpos[i] = pos_cal(xspan, x, xspan_len, xspan_subdiv);
		msh->ypos[i] = pos_cal(yspan, y, yspan_len, yspan_subdiv);
		msh->zpos[i] = pos_cal(zspan, z, zspan_len, zspan_subdiv);
		for(size_t g=0; g<eg_size; ++g){
			msh->chi[g*rt_size+i] = m->chi[g];
			msh->dcoef[g*rt_size+i] = m->dcoef[g];
			msh->sa[g*rt_size+i] = m->sa[g];
			msh->sr[g*rt_size+i] = m->sr[g];
			msh->vsf[g*rt_size+i] = m->vsf[g];
			for(size_t from_g=0; from_g < eg_size; ++from_g)
				msh->ss[g*rt_size+i][from_g] = m->ss[g][from_g];
		}
	}
	
	// ADFs set
	for(size_t g=0; g< eg_size; ++g){
		for(size_t i=0; i< rt_size; ++i){
			msh->adfxl[g*rt_size + i] = 1.0;
			msh->adfxr[g*rt_size + i] = 1.0;
			msh->adfyl[g*rt_size + i] = 1.0;
			msh->adfyr[g*rt_size + i] = 1.0;
			msh->adfzl[g*rt_size + i] = 1.0;
			msh->adfzr[g*rt_size + i] = 1.0;
		}
	}
	for(size_t xspan=0; xspan < xm_span_size; ++xspan){
		for(size_t yspan=0; yspan < ym_span_size; ++yspan){
			for(size_t zspan=0; zspan < zm_span_size; ++zspan){
				if(mtrl_set[xspan][yspan][xspan] < 0)
					continue;
				MBLOCK mblock = sconf_get_mblock(sconf, xspan, yspan, zspan);
				// Set ADF in x direction
				for(size_t j=mblock.start_y; j <= mblock.end_y; ++j)
					for(size_t k=mblock.start_z; k <= mblock.end_z; ++k)
						for(size_t g=0; g<eg_size; ++g){
							////////
						}
				// Set ADF in y direction
				// Set ADF in z direction
			}
		}
	}

	return msh;
}

void msh_free(MSH *msh)
{
	free(msh->mtrl_id);
	free(msh->chi);
	free(msh->dcoef);
	free(msh->sa);
	free(msh->vsf);
	free(msh->dx);
	free(msh->dy);
	free(msh->dz);
	free(msh->xpos);
	free(msh->ypos);
	free(msh->zpos);
	free(msh->adfxl);
	free(msh->adfxr);
	free(msh->adfyl);
	free(msh->adfyr);
	free(msh->adfzl);
	free(msh->adfzr);
	size_t rt_size = msh->rt_size;
	for(size_t i=0; i<rt_size; ++i)
		free(msh->ss[i]);
	free(msh->ss);
	free(msh);
}

inline int msh_get_mtrl_id(const MSH *msh, size_t i, size_t j, size_t k)
{
	return msh->mtrl_id[mapper_get1Didx(msh->mapper, i,j,k)];
}

inline double msh_get_chi(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->chi[g*msh->rt_size + mapper_get1Didx(msh->mapper,i,j,k)];
}

inline double msh_get_dcoef(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->dcoef[g*msh->rt_size + mapper_get1Didx(msh->mapper,i,j,k)];
}

inline double msh_get_sa(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->sa[g*msh->rt_size + mapper_get1Didx(msh->mapper,i,j,k)];
}

inline double msh_get_sr(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->sr[g*msh->rt_size + mapper_get1Didx(msh->mapper,i,j,k)];
}

inline double msh_get_vsf(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->vsf[g*msh->rt_size + mapper_get1Didx(msh->mapper,i,j,k)];
}

inline double msh_get_ss(const MSH *msh, size_t g, size_t from_g, size_t i, size_t j, size_t k)
{
	return msh->ss[g*msh->rt_size + mapper_get1Didx(msh->mapper,i,j,k)][from_g];
}

inline double msh_get_dx(const MSH *msh, size_t i, size_t j, size_t k)
{
	return msh->dx[mapper_get1Didx(msh->mapper,i,j,k)];
}

inline double msh_get_dy(const MSH *msh, size_t i, size_t j, size_t k)
{
	return msh->dy[mapper_get1Didx(msh->mapper,i,j,k)];
}

inline double msh_get_dz(const MSH *msh, size_t i, size_t j, size_t k)
{
	return msh->dz[mapper_get1Didx(msh->mapper,i,j,k)];
}

inline double msh_get_xpos(const MSH *msh, size_t i, size_t j, size_t k)
{
	return msh->xpos[mapper_get1Didx(msh->mapper,i,j,k)];
}

inline double msh_get_ypos(const MSH *msh, size_t i, size_t j, size_t k)
{
	return msh->ypos[mapper_get1Didx(msh->mapper,i,j,k)];
}

inline double msh_get_zpos(const MSH *msh, size_t i, size_t j, size_t k)
{
	return msh->zpos[mapper_get1Didx(msh->mapper,i,j,k)];
}

double msh_get_adfxl(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->adfxl[mapper_get1Didx(msh->mapper,i,j,k)];
}

double msh_get_adfxr(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->adfxr[mapper_get1Didx(msh->mapper,i,j,k)];
}

double msh_get_adfyl(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->adfyl[mapper_get1Didx(msh->mapper,i,j,k)];
}

double msh_get_adfyr(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->adfyr[mapper_get1Didx(msh->mapper,i,j,k)];
}

double msh_get_adfzl(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->adfzl[mapper_get1Didx(msh->mapper,i,j,k)];
}

double msh_get_adfzr(const MSH *msh, size_t g, size_t i, size_t j, size_t k)
{
	return msh->adfzr[mapper_get1Didx(msh->mapper,i,j,k)];
}

void msh_fprintf(const MSH *msh, FILE *stream)
{
	size_t eg_size = msh->eg_size;
	size_t rt_size = msh->rt_size;
	MAPPER *mapper = msh->mapper;
	fprintf(stream, "1DID\t3DXID\t3DYID\t3DZID\tMTRL_ID\tDX\tDY\tDZ\tXPOS\tYPOS\tZPOS\n");
	for(size_t i=0; i< rt_size; ++i){
		XYZ_IDX xyz = mapper->one2three[i];
		size_t r = mapper->three2one[xyz.zi][xyz.yi][xyz.xi];
		fprintf(stream, "%4zd\t%4zd\t%4zd\t%4zd\t%4d\t%4g\t%4g\t%4g\t%4g\t%4g\t%4g\n",
				r, xyz.xi, xyz.yi, xyz.zi, msh->mtrl_id[r], 
				msh->dx[r], msh->dy[r], msh->dz[r], 
				msh->xpos[r], msh->ypos[r], msh->zpos[r]);
	}
	fprintf(stream, "\n");
	fprintf(stream, "1DID\t4DEID\t4DXID\t4DYID\t4DZID\tADFXL\tADFXR\tADFYL\tADFYR\tADFZL\tADFZR\n");
	for(size_t g=0; g< eg_size; ++g){
		for(size_t i=0; i< rt_size; ++i){
			XYZ_IDX xyz = mapper->one2three[i];
			size_t r = mapper->three2one[xyz.zi][xyz.yi][xyz.xi];
			fprintf(stream, "%4zd\t%4zd\t%4zd\t%4zd\t%4zd\t%4g\t%4g\t%4g\t%4g\t%4g\t%4g\n",
				g*rt_size+r, g, xyz.xi, xyz.yi, xyz.zi,
				msh->adfxl[g*rt_size+r], msh->adfxr[g*rt_size+r],
				msh->adfyl[g*rt_size+r], msh->adfyr[g*rt_size+r],
				msh->adfzl[g*rt_size+r], msh->adfzr[g*rt_size+r]);
		}
	}
	fprintf(stream, "\n");
	fprintf(stream, "1DID\t4DEID\t4DXID\t4DYID\t4DZID\tCHI\tDCOEF\tSA\tSR\tVSF\tSS(FROM EG1 ... EGN)\n");
	for(size_t g=0; g< eg_size; ++g){
		for(size_t i=0; i< rt_size; ++i){
			XYZ_IDX xyz = mapper->one2three[i];
			size_t r = mapper->three2one[xyz.zi][xyz.yi][xyz.xi];
			fprintf(stream, "%zd\t%zd\t%zd\t%zd\t%zd\t%g\t%g\t%g\t%g\t%g\t",
					g*rt_size+r, g, xyz.xi, xyz.yi, xyz.zi, msh->chi[g*rt_size+r], 
					msh->dcoef[g*rt_size+r], msh->sa[g*rt_size+r], 
					msh->sr[g*rt_size+r], msh->vsf[g*rt_size+r]);
			for(size_t from_g=0; from_g < eg_size; ++from_g)
				fprintf(stream, "%g\t",msh->ss[g*rt_size+r][from_g]);
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
