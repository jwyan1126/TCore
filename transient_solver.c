#include"transient_solver.h"
#include<math.h>
#include<stdlib.h>

void next_step_cal(FLUX *flux, PCS *pcs, TCONF *tconf, MAPPER *mapper, MESH *mesh, CDAT4 *vsf_lst)
{
	size_t eg_size = mapper->eg_size;
	size_t rt_size = mapper->rt_size;
	CDAT4 *vsf = mesh->vsf;
	FLUX *flux_lst = flux_create(mapper);
	flux_copy(flux_lst, flux);
	CDAT4 *sr_rvs = cdat4_create(mapper);
	sr_rvs_cal(sr_rvs, tconf, mesh);
	MAT *M = mat_create(eg_size*rt_size);
	MAT *S = mat_create(eg_size*rt_size);
	MAT *F = mat_create(eg_size*rt_size);
	EDAT4 *DFDM = edat4_create(mapper);
	EDAT4 *DNOD = edat4_create(mapper);
	cal_DFDM(DFDM, mesh);
	cal_m(M, DFDM, DNOD, mapper, mesh, sr_rvs);
	CDAT4 *chi_bar = cdat4_create(mapper);
	cal_chi_bar(chi_bar, tconf, mesh);
	cal_s(S, mapper, mesh);
	cal_f(F, mapper, mesh, chi_bar);
	mat_adds(M, S, -1.0);
	mat_adds(M, F, -1.0);
	VEC *src_eff = vec_create(eg_size*rt_size);
	src_eff_cal(src_eff,tconf,mesh,vsf_lst,flux_lst,pcs);
	VEC *sol = vec_ref_create(eg_size*rt_size, flux->data);
	(*mat_solver)(sol,M,src_eff,1024);
	vec_ref_free(sol);
	
	pcs_cal(pcs,flux,flux_lst,vsf,vsf_lst,tconf);
	cdat4_free(chi_bar);
	vec_free(src_eff);
	edat4_free(DNOD);
	edat4_free(DFDM);
	mat_free(F);
	mat_free(S);
	mat_free(M);
	cdat4_free(sr_rvs);
	flux_free(flux_lst);
}

void cal_chi_bar(CDAT4 *chi_bar, TCONF *tconf, MESH *mesh)
{
	size_t pcs_size = tconf->pcs_size;
	size_t eg_size = mesh->eg_size;
	size_t rt_size = mesh->rt_size;
	double tau = tconf->tau;
	double *lambdas = tconf->lambdas;
	double *betas = tconf->betas;
	double beta = tconf->beta;
	CDAT4 *chi = mesh->chi;
	MAPPER *mapper = mesh->mapper;
	for(size_t g=0; g<eg_size; ++g)
		for(size_t r=0; r<rt_size; ++r){
			XYZ_IDX xyz = mapper_get3Didx(mapper, r);
			size_t i = xyz.xi; size_t j = xyz.yi; size_t k = xyz.zi;
			double tmp = 0.0;
			for(size_t p=0; p<pcs_size; ++p)
				tmp += betas[p] *(tau - (1.0 - exp(-lambdas[p] * tau))/lambdas[p]) / tau;
			double xg_bar = (1.0-beta)*cdat4_get_val(chi,g,i,j,k) + cdat4_get_val(chi,g,i,j,k)*tmp;
			cdat4_set_val(chi_bar,g,i,j,k,xg_bar);
		}
}

void src_eff_cal(VEC *src_eff, TCONF *tconf, MESH *mesh, CDAT4 *vsf_lst, FLUX *flux_lst, PCS *pcs_lst)
{
	size_t eg_size = mesh->eg_size;
	size_t rt_size = mesh->rt_size;
	size_t pcs_size = tconf->pcs_size;
	CDAT4 *chi = mesh->chi;
	double *lambdas = tconf->lambdas;
	double *betas = tconf->betas;
	double tau = tconf->tau;
	double *nvel = tconf->nvel;
	MAPPER *mapper = mesh->mapper;
	for(size_t g=0; g<eg_size; ++g)
		for(size_t r=0; r<rt_size; ++r){
			XYZ_IDX xyz = mapper_get3Didx(mapper, r);
			size_t i = xyz.xi; size_t j = xyz.yi; size_t k = xyz.zi;
			double sub_term1 = flux_lst->data[g*rt_size+r]/(nvel[g]*tau);
			double sub_term2 = 0.0;
			for(size_t p=0; p<pcs_size; ++p)
				sub_term2 += lambdas[p] * pcs_lst->data[p][r] * exp(-lambdas[p]*tau);
			sub_term2 *= cdat4_get_val(chi,g,i,j,k);

			double sub_term3 = 0.0;
			for(size_t p=0; p<pcs_size; ++p){
				double tmp = 0.0;
				for(size_t from_g=0; from_g<eg_size; ++from_g)
					tmp += cdat4_get_val(vsf_lst,from_g,i,j,k) * flux_lst->data[from_g*rt_size+r];
				double F1 = betas[p] * (tau - (1.0 - exp(-lambdas[p] * tau))/ lambdas[p]) / (tau * lambdas[p]);
				double F0 = betas[p] * (1.0 - exp(-lambdas[p] * tau)) / lambdas[p] - F1;
				tmp *= lambdas[p];
				tmp *= F0;
				sub_term3 += tmp;
			}
			sub_term3 *= cdat4_get_val(chi,g,i,j,k);
			src_eff->vals[g*rt_size+r] = sub_term1 + sub_term2 + sub_term3;
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
	int ***cchecker = sr_rvs->cchecker;
	for(size_t k=0; k<zm_size; ++k)
		for(size_t j=0; j<ym_size; ++j)
			for(size_t i=0; i<xm_size; ++i){
				if(cchecker[k][j][i] & 0b00000001) continue;
				for(size_t g=0; g<eg_size; ++g){
					double val = cdat4_get_val(sr,g,i,j,k) + 1.0/(nvel[g]*tau);
					cdat4_set_val(sr_rvs,g,i,j,k, val);
				}
			}
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
			for(size_t from_g=0; from_g<eg_size; ++from_g){
				tmp1 += cdat4_get_val(vsf_lst,from_g,i,j,k) * flux_lst->data[from_g*rt_size+r];
				tmp2 += cdat4_get_val(vsf,from_g,i,j,k) * flux->data[from_g*rt_size+r];
			}
			pcs->data[p][r] = pcs->data[p][r]*exp(-lambdas[p]*tau) + F0*tmp1 + F1*tmp2;
		}
	}
}
