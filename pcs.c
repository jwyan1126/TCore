#include"pcs.h"
#include<stdlib.h>

PCS *pcs_create(size_t pcs_size, MAPPER *mapper)
{
	PCS *pcs = malloc(sizeof(PCS));
	pcs->pcs_size = pcs_size;
	pcs->xm_size = mapper->xm_size;
	pcs->ym_size = mapper->ym_size;
	pcs->zm_size = mapper->zm_size;
	size_t rt_size = mapper->rt_size;
	pcs->rt_size = rt_size;
	pcs->mapper = mapper;
	pcs->data = malloc(pcs_size * sizeof(double *));
	for(size_t p=0; p<pcs_size; ++p)
		pcs->data[p] = calloc(rt_size, sizeof(double));
	return pcs;
}

void pcs_free(PCS *pcs)
{
	size_t pcs_size = pcs->pcs_size;
	for(size_t p=0; p<pcs_size; ++p)
		free(pcs->data[p]);
	free(pcs->data);
	free(pcs);
}

double pcs_get_val(const PCS *pcs, size_t p, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	size_t pcs_size = pcs->pcs_size;
	size_t xm_size = pcs->xm_size;
	size_t ym_size = pcs->ym_size;
	size_t zm_size = pcs->zm_size;
	int ***cchecker = pcs->mapper->cchecker;
	if(p >= pcs_size || i >= xm_size || j >= ym_size || k >= zm_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	if(cchecker[k][j][i] & 0b00000001){
		fprintf(stderr, "No mtrl filled in.\n");
		exit(-1);
	}
	#endif
	size_t idx = mapper_get1Didx(pcs->mapper, i, j, k);
	return pcs->data[p][idx];
}

void pcs_copy(PCS *tar_pcs, const PCS *src_pcs)
{
	#ifdef DEBUG
	if(tar_pcs->pcs_size != src_pcs->pcs_size ||
	   tar_pcs->xm_size != src_pcs->xm_size ||
	   tar_pcs->ym_size != src_pcs->ym_size ||
	   tar_pcs->zm_size != src_pcs->zm_size ||	
	   tar_pcs->rt_size != src_pcs->rt_size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	tar_pcs->pcs_size = src_pcs->pcs_size;
	tar_pcs->xm_size = src_pcs->xm_size;
	tar_pcs->ym_size = src_pcs->ym_size;
	tar_pcs->zm_size = src_pcs->zm_size;
	tar_pcs->rt_size = src_pcs->rt_size;
	tar_pcs->mapper = src_pcs->mapper;
	for(size_t p=0; p<tar_pcs->pcs_size; ++p)
		for(size_t i=0; i<tar_pcs->rt_size; ++i)
			tar_pcs->data[p][i] = src_pcs->data[p][i];
}

void pcs_init(PCS *pcs, FLUX *flux, MESH *mesh, TCONF *tconf)
{
	size_t pcs_size = pcs->pcs_size;
	size_t eg_size = mesh->eg_size;
	size_t rt_size = pcs->rt_size;
	MAPPER *mapper = pcs->mapper;
	double *betas = tconf->betas;
	double *lambdas = tconf->lambdas;
	for(size_t p=0; p< pcs_size; ++p){
		for(size_t r=0; r< rt_size; ++r){
			XYZ_IDX xyz = mapper_get3Didx(mapper, r);
			size_t i = xyz.xi; size_t j = xyz.yi; size_t k = xyz.zi;
			double sum = 0.0;
			for(size_t g=0; g< eg_size; ++g)
				sum += cdat4_get_val(mesh->vsf,g,i,j,k) * flux->data[g*rt_size+r];
			sum *= betas[p];
			pcs->data[p][r] = sum / lambdas[p];
		}
	}
}

void pcs_fprintf(const PCS *pcs, FILE *stream)
{
	size_t pcs_size = pcs->pcs_size;
	size_t xm_size = pcs->xm_size;
	size_t ym_size = pcs->ym_size;
	size_t zm_size = pcs->zm_size;
	MAPPER *mapper = pcs->mapper;
	int ***cchecker = mapper->cchecker;
	for(size_t p=0; p<pcs_size; ++p){
		fprintf(stream, "p=%zd\n", p);
		for(size_t k=0; k<zm_size; ++k){
			fprintf(stream, "k=%zd\n", k);
			for(size_t j=0; j<ym_size; ++j){
				for(size_t i=0; i<xm_size; ++i){
					size_t idx = mapper_get1Didx(mapper, i, j, k);
					double val = (cchecker[k][j][i] & 0b00000001) ? 0.0 : pcs->data[p][idx];
					fprintf(stream, "%g\t", val);
				}
				fprintf(stream, "\n");
			}
			fprintf(stream, "\n");
		}
		fprintf(stream, "\n");
	}
}
