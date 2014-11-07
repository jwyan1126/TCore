#include"tsol.h"
#include<stdlib.h>


TSOL *tsol_create(const TCONF *tconf, const MESH *mesh, const SSOL *ssol)
{
	MAPPER *mapper = ssol->flux->mapper;
	size_t pcs_size = tconf->pcs_size;
	size_t eg_size = mapper->eg_size;
	size_t rt_size = mapper->rt_size;
	TSOL *tsol = malloc(sizeof(TSOL));
	tsol->eg_size = eg_size;
	tsol->rt_size = rt_size;
	tsol->xm_size = mapper->xm_size;
	tsol->ym_size = mapper->ym_size;
	tsol->zm_size = mapper->zm_size;
	tsol->pcs_size = pcs_size;
	tsol->time = 0.0;
	tsol->keff = ssol->keff;
	tsol->flux = flux_create(mapper);
	flux_copy(tsol->flux, ssol->flux);
	tsol->pcs = pcs_create(pcs_size, mapper);
	// cal init pcs
	for(size_t p=0; p<pcs_size; ++p)
		for(size_t r=0; r<rt_size; ++r){
			XYZ_IDX xyz = mapper_get3Didx(mapper, r);
			double sum = 0.0;
			for(size_t g=0; g<eg_size; ++g)
				sum += cdat4_get_val(mesh->vsf,g,xyz.xi,xyz.yi,xyz.zi) * tsol->flux->data[g*rt_size+r];
			sum *= tconf->betas[p];
			tsol->pcs->data[p][r] = sum / tconf->lambdas[p];
		}
	return tsol;
}

void tsol_free(TSOL *tsol)
{
	flux_free(tsol->flux);
	pcs_free(tsol->pcs);
	free(tsol);
}
