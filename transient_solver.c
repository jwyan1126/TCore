#include"transient_solver.h"
#include<math.h>
#include<stdlib.h>

void next_step_cal(TSOL *tsol, TCONF *tconf, MAPPER *mapper, MESH *mesh)
{
	tsol->time += tconf->tau;
	size_t eg_size = mapper->eg_size;
	size_t rt_size = mapper->rt_size;
	FLUX *flux = tsol->flux;
	PCS *pcs = tsol->pcs;
	CDAT4 *vsf = mesh->vsf;
	static FLUX *flux_lst;
	static CDAT4 *vsf_lst;
	static PCS *pcs_lst;
	if(flux_lst == NULL) flux_lst = flux_create(mapper);
	if(vsf_lst == NULL) vsf_lst = cdat4_create(mapper);
	if(pcs_lst == NULL) pcs_lst = pcs_create(tsol->pcs_size, mapper);
	flux_copy(flux_lst, flux);
	cdat4_copy(vsf_lst, vsf);
	pcs_copy(pcs_lst, pcs);
	static CDAT4 *sr_rvs;
	if(sr_rvs == NULL) sr_rvs = cdat4_create(mapper);
	sr_rvs_cal(sr_rvs, tconf, mesh);
	static MAT *M;
	if(M == NULL) M = mat_create(eg_size*rt_size);
	static EDAT4 *DFDM;
	if(DFDM == NULL) DFDM = edat4_create(mapper);
	static EDAT4 *DNOD;
	if(DNOD == NULL) DNOD = edat4_create(mapper);
	cal_DFDM(DFDM, mesh);
	cal_m(M, DFDM, DNOD, mapper, mesh, sr_rvs);
	static VEC *src_rvs;
	if(src_rvs == NULL) src_rvs = vec_create(eg_size*rt_size);
	src_rvs_cal(src_rvs,tconf,mesh,vsf_lst,flux,flux_lst,pcs,pcs_lst);
	VEC *sol = vec_ref_create(eg_size*rt_size, flux->data);
	(*mat_solver)(sol,M,src_rvs,256);
	vec_ref_free(sol);
	pcs_cal(pcs,flux,flux_lst,vsf,vsf_lst,tconf);
}

void src_rvs_cal(VEC *src, TCONF *tconf, MESH *mesh, CDAT4 *vsf_lst, FLUX *flux, FLUX *flux_lst, PCS *pcs, PCS *pcs_lst)
{
	size_t eg_size = mesh->eg_size;
	size_t rt_size = mesh->rt_size;
	size_t pcs_size = tconf->pcs_size;
	CDAT5 *ss = mesh->ss;
	CDAT4 *vsf = mesh->vsf;
	CDAT4 *chi = mesh->chi;
	double *lambdas = tconf->lambdas;
	double *betas = tconf->betas;
	double beta = tconf->beta;
	double tau = tconf->tau;
	double *nvel = tconf->nvel;
	MAPPER *mapper = mesh->mapper;
	for(size_t g=0; g<eg_size; ++g)
		for(size_t r=0; r<rt_size; ++r){
			XYZ_IDX xyz = mapper_get3Didx(mapper, r);
			size_t i = xyz.xi; size_t j = xyz.yi; size_t k = xyz.zi;
			double sca_term, fis_term, eff_term;
			// cal. sca_term
			sca_term = 0.0;
			for(size_t from_g=0; from_g<eg_size; ++from_g)
				sca_term += cdat5_get_val(ss,g,from_g,i,j,k) * flux->data[from_g*rt_size+r];
			// cal. fis_term
			double tmp = 0.0;
			for(size_t p=0; p<pcs_size; ++p)
				tmp += lambdas[p] * betas[p] *(tau - (1.0 - exp(-lambdas[p] * tau))/lambdas[p]) / (tau*lambdas[p]);
			double xg_bar = (1.0-beta)*cdat4_get_val(chi,g,i,j,k) + cdat4_get_val(chi,g,i,j,k)*tmp;
			fis_term = 0.0;
			for(size_t from_g; from_g<eg_size; ++from_g)
				fis_term += cdat4_get_val(vsf,from_g,i,j,k) * flux->data[from_g*rt_size+r];
			fis_term *= xg_bar;
			// cal. eff_term
			double sub_term1 = flux_lst->data[g*rt_size+r]/(nvel[g]*tau);
			double sub_term2 = 0.0;
			for(size_t p=0; p<pcs_size; ++p)
				sub_term2 += lambdas[p] * pcs_lst->data[p][r] * exp(-lambdas[p]*tau);
			sub_term2 *= cdat4_get_val(chi,g,i,j,k);

			double sub_term3 = 0.0;
			for(size_t p=0; p<pcs_size; ++p){
				tmp = 0.0;
				for(size_t from_g=0; from_g<eg_size; ++from_g)
					tmp += cdat4_get_val(vsf_lst,from_g,i,j,k) * flux_lst->data[from_g*rt_size+r];
				double F1 = betas[p] * (tau - (1.0 - exp(-lambdas[p] * tau))/ lambdas[p]) / (tau * lambdas[p]);
				double F0 = betas[p] * (1.0 - exp(-lambdas[p] * tau)) / lambdas[p] - F1;
				tmp *= lambdas[p];
				tmp *= F0;
				sub_term3 += tmp;
			}
			sub_term3 *= cdat4_get_val(chi,g,i,j,k);
			eff_term = sub_term1 + sub_term2 + sub_term3;
			src->vals[g*rt_size+r] = sca_term + fis_term + eff_term;
		}
}

void sr_rvs_cal(CDAT4 *sr_rvs, TCONF *tconf, MESH *mesh)
{
	CDAT4 *sr = mesh->sr;
	double *nvel = tconf->nvel;
	double tau = tconf->tau;
	size_t eg_size = mesh->eg_size;
	size_t xm_size = mesh->xm_mesh_size;
	size_t ym_size = mesh->ym_mesh_size;
	size_t zm_size = mesh->zm_mesh_size;
	for(size_t k=0; k<zm_size; ++k)
		for(size_t j=0; j<ym_size; ++j)
			for(size_t i=0; i<xm_size; ++i)
				for(size_t g=0; g<eg_size; ++g)
					sr_rvs->data[k][j][i][g] = sr->data[k][j][i][g] + 1.0/(nvel[g] * tau);
}

void pcs_cal(PCS *pcs, FLUX *flux, FLUX *flux_lst, CDAT4 *vsf, CDAT4 *vsf_lst, TCONF *tconf)
{
	MAPPER *mapper = flux->mapper;
	size_t eg_size = mapper->eg_size;
	size_t pcs_size = pcs->pcs_size;
	size_t rt_size = mapper->rt_size;
	double *lambdas = tconf->lambdas;
	double *betas = tconf->betas;
	double tau = tconf->tau;
	for(size_t p=0; p<pcs_size; ++p){
		double F1 = betas[p] * (tau - (1.0 - exp(-lambdas[p] * tau))/ lambdas[p]) / (tau * lambdas[p]);
		double F0 = betas[p] * (1.0 - exp(-lambdas[p] * tau)) / lambdas[p] - F1;
		for(size_t r=0; r<rt_size; ++r){
			XYZ_IDX xyz = mapper_get3Didx(mapper, r);
			size_t i = xyz.xi; size_t j = xyz.yi; size_t k = xyz.zi;
			double tmp1 = 0.0;
			double tmp2 = 0.0;
			for(size_t from_g; from_g<eg_size; ++from_g){
				tmp1 += cdat4_get_val(vsf_lst,from_g,i,j,k) * flux_lst->data[from_g*rt_size+r];
				tmp2 += cdat4_get_val(vsf,from_g,i,j,k) * flux->data[from_g*rt_size+r];
			}
			pcs->data[p][r] = pcs->data[p][r]*exp(-lambdas[p]*tau) + F0*tmp1 + F1*tmp2;
		}
	}
}
