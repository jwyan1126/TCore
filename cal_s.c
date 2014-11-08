#include"steady_solver.h"

void cal_s(MAT *S, const MAPPER *mapper, const MESH *mesh)
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
