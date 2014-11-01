#include"leak.h"

LEAK *leak_create(MAPPER *mapper)
{
	LEAK *leak = malloc(sizeof(LEAK));
	leak->eg_size = mapper->eg_size;
	leak->xm_size = mapper->xm_size;
	leak->ym_size = mapper->ym_size;
	leak->zm_size = mapper->zm_size;
	leak->rt_size = mapper->rt_size;
	leak->lx0 = cdat4_create(mapper);
	leak->lx1 = cdat4_create(mapper);
	leak->lx2 = cdat4_create(mapper);
	leak->ly0 = cdat4_create(mapper);
	leak->ly1 = cdat4_create(mapper);
	leak->ly2 = cdat4_create(mapper);
	leak->lz0 = cdat4_create(mapper);
	leak->lz1 = cdat4_create(mapper);
	leak->lz2 = cdat4_create(mapper);
	leak->cchecker = mapper->cchecker;
	return leak;
}

void leak_free(LEAK *leak)
{
	cdat4_free(lx0);
	cdat4_free(lx1);
	cdat4_free(lx2);
	cdat4_free(ly0);
	cdat4_free(ly1);
	cdat4_free(ly2);
	cdat4_free(lz0);
	cdat4_free(lz1);
	cdat4_free(lz2);
	free(leak);
}

void cal_leakage(LEAK *leak, const MESH *mesh, const EDAT4 *jn)
{
	size_t eg_size = jn->eg_size;
	size_t xm_size = jn->xm_size;
	size_t ym_size = jn->ym_size;
	size_t zm_size = jn->zm_size;
	for(size_t k=0; k<zm_size; ++k)
		for(size_t j=0; j<ym_size; ++j)
			for(size_t i=0; i<xm_size; ++i){
				double dx = cdat3_get_val(mesh->dx, i, j, k);
				double dy = cdat3_get_val(mesh->dy, i, j, k);
				double dz = cdat3_get_val(mesh->dz, i, j, k);
				for(size_t g=0; g<eg_size; ++g){
					double Jgxl = edat4_get_xlval(jn, g, i, j, k);
					double Jgxr = edat4_get_xrval(jn, g, i, j, k);
					double Jgyl = edat4_get_ylval(jn, g, i, j, k);
					double Jgyr = edat4_get_yrval(jn, g, i, j, k);
					double Jgzl = edat4_get_zlval(jn, g, i, j, k);
					double Jgzr = edat4_get_zrval(jn, g, i, j, k);
					double lgx0 = (Jgyr - Jgyl) / dy + (Jgzr - Jgzl) / dz;
					double lgy0 = (Jgzr - Jgzl) / dz + (Jgxr - Jgxl) / dx;
					double lgz0 = (Jgxr - Jgxl) / dx + (Jgyr - Jgyl) / dy;
					cdat4_set_val(leak->lx0, g, i, j, k, lgx0);
					cdat4_set_val(leak->ly0, g, i, j, k, lgy0);
					cdat4_set_val(leak->lz0, g, i, j, k, lgz0);
				}
			}
}
