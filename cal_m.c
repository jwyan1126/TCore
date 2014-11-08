#include"steady_solver.h"

void cal_m(MAT *M, EDAT4 *DFDM, EDAT4 *DNOD, const MAPPER *mapper, const MESH *mesh, CDAT4 *sr_rvs)
{
	size_t eg_size = mesh->eg_size;
	size_t rt_size = mapper->rt_size;
	CDAT3 *dx = mesh->dx;
	CDAT3 *dy = mesh->dy;
	CDAT3 *dz = mesh->dz;
	CDAT4 *sr;
	if(sr_rvs == NULL) sr = mesh->sr;
	else sr = sr_rvs;
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
			mat_h += (edat4_get_xrval(DFDM,g,i,j,k) - edat4_get_xrval(DNOD,g,i,j,k)) 
				* cdat4_get_val(adfxr,g,i,j,k) / cdat3_get_val(dx,i,j,k);
			mat_h += (edat4_get_xlval(DFDM,g,i,j,k) + edat4_get_xlval(DNOD,g,i,j,k)) 
				* cdat4_get_val(adfxl,g,i,j,k) / cdat3_get_val(dx,i,j,k);
			mat_h += (edat4_get_yrval(DFDM,g,i,j,k) - edat4_get_yrval(DNOD,g,i,j,k)) 
				* cdat4_get_val(adfyr,g,i,j,k) / cdat3_get_val(dy,i,j,k);
			mat_h += (edat4_get_ylval(DFDM,g,i,j,k) + edat4_get_ylval(DNOD,g,i,j,k)) 
				* cdat4_get_val(adfyl,g,i,j,k) / cdat3_get_val(dy,i,j,k);
			mat_h += (edat4_get_zrval(DFDM,g,i,j,k) - edat4_get_zrval(DNOD,g,i,j,k)) 
				* cdat4_get_val(adfzr,g,i,j,k) / cdat3_get_val(dz,i,j,k);
			mat_h += (edat4_get_zlval(DFDM,g,i,j,k) + edat4_get_zlval(DNOD,g,i,j,k)) 
				* cdat4_get_val(adfzl,g,i,j,k) / cdat3_get_val(dz,i,j,k);
			mat_set(M, p, p, mat_h);
			
			// XL  
			if(!(checker[k][j][i] & 0b00000100)){
				double mat_b = -(edat4_get_xlval(DFDM,g,i,j,k) - edat4_get_xlval(DNOD,g,i,j,k)) 
						* cdat4_get_val(adfxr,g,i-1,j,k) / cdat3_get_val(dx,i,j,k);
				size_t pxl = g*rt_size + mapper_get1Didx(mapper,i-1,j,k);
				mat_set(M, p, pxl, mat_b);
			}
			// XR
			if(!(checker[k][j][i] & 0b00001000)){
				double mat_a = -(edat4_get_xrval(DFDM,g,i,j,k) + edat4_get_xrval(DNOD,g,i,j,k)) 
						* cdat4_get_val(adfxl,g,i+1,j,k) / cdat3_get_val(dx,i,j,k);
				size_t pxr = g*rt_size + mapper_get1Didx(mapper,i+1,j,k);
				mat_set(M, p, pxr, mat_a);
			}
			// YL
			if(!(checker[k][j][i] & 0b00010000)){
				double mat_d = -(edat4_get_ylval(DFDM,g,i,j,k) - edat4_get_ylval(DNOD,g,i,j,k)) 
						* cdat4_get_val(adfyr,g,i,j-1,k) / cdat3_get_val(dy,i,j,k);
				size_t pyl = g*rt_size + mapper_get1Didx(mapper,i,j-1,k);
				mat_set(M, p, pyl, mat_d);
			}
			// YR
			if(!(checker[k][j][i] & 0b00100000)){
				double mat_c = -(edat4_get_yrval(DFDM,g,i,j,k) + edat4_get_yrval(DNOD,g,i,j,k)) 
						* cdat4_get_val(adfyl,g,i,j+1,k) / cdat3_get_val(dy,i,j,k);
				size_t pyr = g*rt_size + mapper_get1Didx(mapper,i,j+1,k);
				mat_set(M, p, pyr, mat_c);
			}
			// ZL
			if(!(checker[k][j][i] & 0b01000000)){
				double mat_f = -(edat4_get_zlval(DFDM,g,i,j,k) - edat4_get_zlval(DNOD,g,i,j,k)) 
						* cdat4_get_val(adfzr,g,i,j,k-1) / cdat3_get_val(dz,i,j,k);
				size_t pzl = g*rt_size + mapper_get1Didx(mapper,i,j,k-1);
				mat_set(M, p, pzl, mat_f);
			}
			// ZR
			if(!(checker[k][j][i] & 0b10000000)){
				double mat_e = -(edat4_get_zrval(DFDM,g,i,j,k) + edat4_get_zrval(DNOD,g,i,j,k)) 
						* cdat4_get_val(adfzl,g,i,j,k+1) / cdat3_get_val(dz,i,j,k);
				size_t pzr = g*rt_size + mapper_get1Didx(mapper,i,j,k+1);
				mat_set(M, p, pzr, mat_e);
			}
		}
	}
}
