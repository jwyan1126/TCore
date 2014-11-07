#include"steady_solver.h"

void cal_f(MAT *F, const SCONF *sconf, const MAPPER *mapper, const MESH *mesh)
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
