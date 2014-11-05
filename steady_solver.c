#include"steady_solver.h"
#include"ssol.h"
#include"algebra/mat.h"
#include"algebra/vec.h"
#include"algebra/eigen.h"
#include<stdlib.h>
#include"nodal/leak.h"
#include"nodal/jcur.h"
#include"nodal/tnsol.h"
#include"nodal/SANM/sanm.h"

void steady_solver(SSOL *ssol, SCONF *sconf, MAPPER *mapper, const MESH *mesh)
{
	size_t eg_size = sconf->eg_size;
	size_t rt_size = sconf->rt_mesh_size;
	MAT *M = mat_create(eg_size * rt_size);
	MAT *S = mat_create(eg_size * rt_size);
	MAT *F = mat_create(eg_size * rt_size);
	EDAT4 *DFDM = edat4_create(mapper);
	EDAT4 *DNOD = edat4_create(mapper);
	EDAT4 *Jn = edat4_create(mapper);
	LEAK *leak = leak_create(mapper);
	VEC *phi = vec_ref_create(eg_size * rt_size, ssol->flux->data);
	vec_normalize(phi);
	VEC *phi_last = vec_create(eg_size * rt_size);
	
	cal_DFDM(DFDM, mesh);
	cal_S(S, sconf,mapper, mesh);
	cal_F(F, sconf, mapper, mesh);
	double res = 1.0;
	int counter = 0;
	//while(res > 1e-6){
	for(int ggg =0; ggg < 10000; ++ggg){
		cal_M(M, DFDM, DNOD, mapper, mesh);
		mat_adds(M, S, -1.0);
		double k;
		vec_copy(phi_last, phi);
		gspow_iter(&k, phi, M, F, 0.0, 4096);
		vec_abs(phi);
		vec_normalize(phi);
		res = vec_res_2norm(phi_last, phi);
		counter++;
		ssol->keff = 1.0 / k;
		printf("nonlinear iter = %d\tres = %g\tkeff = %g\n", counter, res, ssol->keff);
		cal_jcur(Jn, mesh, leak, ssol);
		cal_leakage(leak, mesh, Jn);
		cal_DNOD(DNOD, mesh, DFDM, Jn, ssol);
	}
	printf("keff=%g\n",ssol->keff);
	vec_free(phi_last);
	leak_free(leak);
	edat4_free(Jn);
	edat4_free(DNOD);
	edat4_free(DFDM);
	vec_ref_free(phi);
	mat_free(F);
	mat_free(S);
	mat_free(M);
}

void cal_M(MAT *M, EDAT4 *DFDM, EDAT4 *DNOD, const MAPPER *mapper, const MESH *mesh)
{
	size_t eg_size = mesh->eg_size;
	size_t rt_size = mapper->rt_size;
	CDAT3 *dx = mesh->dx;
	CDAT3 *dy = mesh->dy;
	CDAT3 *dz = mesh->dz;
	CDAT4 *sr = mesh->sr;
	CDAT4 *adfxl = mesh->adfxl;
	CDAT4 *adfxr = mesh->adfxr;
	CDAT4 *adfyl = mesh->adfyl;
	CDAT4 *adfyr = mesh->adfyr;
	CDAT4 *adfzl = mesh->adfzl;
	CDAT4 *adfzr = mesh->adfzr;
	int ***checker = mesh->cchecker;
	mat_set_zeros(M);
	for(size_t g=0; g<eg_size; ++g){
		for(size_t idx=0; idx<rt_size; ++idx){
			XYZ_IDX xyz = mapper_get3Didx(mapper, idx);
			size_t i = xyz.xi;
			size_t j = xyz.yi;
			size_t k = xyz.zi;
			size_t p = g*rt_size + idx;
			double mat_h = 0.0;
			mat_h += cdat4_get_val(sr, g, i, j, k);
			mat_h += (edat4_get_xrval(DFDM,g,i,j,k) - edat4_get_xrval(DNOD,g,i,j,k)) * cdat4_get_val(adfxr,g,i,j,k) / cdat3_get_val(dx,i,j,k);
			mat_h += (edat4_get_xlval(DFDM,g,i,j,k) + edat4_get_xlval(DNOD,g,i,j,k)) * cdat4_get_val(adfxl,g,i,j,k) / cdat3_get_val(dx,i,j,k);
			mat_h += (edat4_get_yrval(DFDM,g,i,j,k) - edat4_get_yrval(DNOD,g,i,j,k)) * cdat4_get_val(adfyr,g,i,j,k) / cdat3_get_val(dy,i,j,k);
			mat_h += (edat4_get_ylval(DFDM,g,i,j,k) + edat4_get_ylval(DNOD,g,i,j,k)) * cdat4_get_val(adfyl,g,i,j,k) / cdat3_get_val(dy,i,j,k);
			mat_h += (edat4_get_zrval(DFDM,g,i,j,k) - edat4_get_zrval(DNOD,g,i,j,k)) * cdat4_get_val(adfzr,g,i,j,k) / cdat3_get_val(dz,i,j,k);
			mat_h += (edat4_get_zlval(DFDM,g,i,j,k) + edat4_get_zlval(DNOD,g,i,j,k)) * cdat4_get_val(adfzl,g,i,j,k) / cdat3_get_val(dz,i,j,k);
			mat_set(M, p, p, mat_h);
			
			// XL  
			if(!(checker[k][j][i] & 0b00000100)){
				double mat_b = -(edat4_get_xlval(DFDM,g,i,j,k) - edat4_get_xlval(DNOD,g,i,j,k)) * cdat4_get_val(adfxr,g,i-1,j,k) / cdat3_get_val(dx,i,j,k);
				size_t pxl = g*rt_size + mapper_get1Didx(mapper,i-1,j,k);
				mat_set(M, p, pxl, mat_b);
			}
			// XR
			if(!(checker[k][j][i] & 0b00001000)){
				double mat_a = -(edat4_get_xrval(DFDM,g,i,j,k) + edat4_get_xrval(DNOD,g,i,j,k)) * cdat4_get_val(adfxl,g,i+1,j,k) / cdat3_get_val(dx,i,j,k);
				size_t pxr = g*rt_size + mapper_get1Didx(mapper,i+1,j,k);
				mat_set(M, p, pxr, mat_a);
			}
			// YL
			if(!(checker[k][j][i] & 0b00010000)){
				double mat_d = -(edat4_get_ylval(DFDM,g,i,j,k) - edat4_get_ylval(DNOD,g,i,j,k)) * cdat4_get_val(adfyr,g,i,j-1,k) / cdat3_get_val(dy,i,j,k);
				size_t pyl = g*rt_size + mapper_get1Didx(mapper,i,j-1,k);
				mat_set(M, p, pyl, mat_d);
			}
			// YR
			if(!(checker[k][j][i] & 0b00100000)){
				double mat_c = -(edat4_get_yrval(DFDM,g,i,j,k) + edat4_get_yrval(DNOD,g,i,j,k)) * cdat4_get_val(adfyl,g,i,j+1,k) / cdat3_get_val(dy,i,j,k);
				size_t pyr = g*rt_size + mapper_get1Didx(mapper,i,j+1,k);
				mat_set(M, p, pyr, mat_c);
			}
			// ZL
			if(!(checker[k][j][i] & 0b01000000)){
				double mat_f = -(edat4_get_zlval(DFDM,g,i,j,k) - edat4_get_zlval(DNOD,g,i,j,k)) * cdat4_get_val(adfzr,g,i,j,k-1) / cdat3_get_val(dz,i,j,k);
				size_t pzl = g*rt_size + mapper_get1Didx(mapper,i,j,k-1);
				mat_set(M, p, pzl, mat_f);
			}
			// ZR
			if(!(checker[k][j][i] & 0b10000000)){
				double mat_e = -(edat4_get_zrval(DFDM,g,i,j,k) + edat4_get_zrval(DNOD,g,i,j,k)) * cdat4_get_val(adfzl,g,i,j,k+1) / cdat3_get_val(dz,i,j,k);
				size_t pzr = g*rt_size + mapper_get1Didx(mapper,i,j,k+1);
				mat_set(M, p, pzr, mat_e);
			}
		}
	}
}

void cal_S(MAT *S, const SCONF *sconf, const MAPPER *mapper, const MESH *mesh)
{
	CDAT5 *ss = mesh->ss;
	mat_set_zeros(S);
	size_t eg_size = mesh->eg_size;
	size_t rt_size = mapper->rt_size;
	for(size_t g=0; g<eg_size; ++g){
		for(size_t idx=0; idx<rt_size; ++idx){
			for(size_t from_g=0; from_g<eg_size; ++from_g){
				if(from_g == g) continue;
				size_t p = g*rt_size + idx;
				size_t from_p = from_g*rt_size + idx;
				XYZ_IDX xyz = mapper_get3Didx(mapper, idx);
				mat_set(S, p, from_p, cdat5_get_val(ss, g, from_g, xyz.xi, xyz.yi, xyz.zi));
			}
		}
	}
}

void cal_F(MAT *F, const SCONF *sconf, const MAPPER *mapper, const MESH *mesh)
{
	CDAT4 *vsf = mesh->vsf;
	CDAT4 *chi = mesh->chi;
	mat_set_zeros(F);
	size_t eg_size = mesh->eg_size;
	size_t rt_size = mapper->rt_size;
	for(size_t g=0; g<eg_size; ++g){
		for(size_t idx=0; idx<rt_size; ++idx){
			for(size_t from_g=0; from_g<eg_size; ++from_g){
				size_t p = g*rt_size + idx;
				size_t from_p = from_g*rt_size + idx;
				XYZ_IDX xyz = mapper_get3Didx(mapper, idx);
				double val = cdat4_get_val(chi, g, xyz.xi, xyz.yi, xyz.zi) * cdat4_get_val(vsf, from_g, xyz.xi, xyz.yi, xyz.zi);
				mat_set(F, p, from_p, val);
			}
		}
	}
}
