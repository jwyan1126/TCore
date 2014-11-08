#include"flux.h"
#include"stdlib.h"

FLUX *flux_create(MAPPER *mapper)
{
	FLUX *flux = malloc(sizeof(FLUX));
	size_t eg_size = mapper->eg_size;
	size_t rt_size = mapper->rt_size;
	flux->eg_size = eg_size;
	flux->rt_size = rt_size;
	flux->xm_size = mapper->xm_size;
	flux->ym_size = mapper->ym_size;
	flux->zm_size = mapper->zm_size;
	flux->data = calloc(rt_size * eg_size, sizeof(double));
	for(size_t i=0; i<eg_size*rt_size; ++i)
		flux->data[i] = 1.0;
	flux->mapper = mapper;
	return flux;
}

void flux_free(FLUX *flux)
{
	free(flux->data);
	free(flux);
}

void flux_copy(FLUX *tar_flux, const FLUX *src_flux)
{
	#ifdef DEBUG
	if(tar_flux->eg_size != src_flux->eg_size ||
	   tar_flux->xm_size != src_flux->xm_size ||
	   tar_flux->ym_size != src_flux->ym_size ||
	   tar_flux->zm_size != src_flux->zm_size ||
	   tar_flux->rt_size != src_flux->rt_size){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	size_t eg_size = src_flux->eg_size;
	size_t rt_size = src_flux->rt_size;
	tar_flux->eg_size = eg_size;
	tar_flux->rt_size = rt_size;
	tar_flux->xm_size = src_flux->xm_size;
	tar_flux->ym_size = src_flux->ym_size;
	tar_flux->zm_size = src_flux->zm_size;
	for(size_t i=0; i< eg_size * rt_size; ++i)
		tar_flux->data[i] = src_flux->data[i];
	tar_flux->mapper = src_flux->mapper;
}

void flux_fprintf(const FLUX *flux, FILE *stream)
{
	MAPPER *mapper = flux->mapper;
	size_t eg_size = mapper->eg_size;
	size_t rt_size = mapper->rt_size;
	size_t xm_size = mapper->xm_size;
	size_t ym_size = mapper->ym_size;
	size_t zm_size = mapper->zm_size;
	int ***cchecker = mapper->cchecker;
	for(size_t g=0; g<eg_size; ++g){
		fprintf(stream, "g=%zd\n", g);
		for(size_t k=0; k<zm_size; ++k){
			fprintf(stream, "k=%zd\n", k);
			for(size_t j=0; j<ym_size; ++j){
				for(size_t i=0; i<xm_size; ++i){
					size_t idx = g*rt_size + mapper_get1Didx(mapper, i, j, k);
					double val = (cchecker[k][j][i] & 0b00000001) ? 0.0 : flux->data[idx];
					fprintf(stream, "%g\t", val);
				}
				fprintf(stream, "\n");
			}
			fprintf(stream, "\n");
		}
		fprintf(stream, "\n");
	}
}

void flux_normalize(FLUX *flux)
{
	size_t eg_size = flux->eg_size;
	size_t rt_size = flux->rt_size;
	double s = 0.0;
	for(size_t i=0; i<rt_size * eg_size; ++i)
		s += flux->data[i];
	s /= (rt_size * eg_size);
	for(size_t i=0; i<rt_size * eg_size; ++i)
		flux->data[i] /= s;
}

double flux_get_val(FLUX *flux, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	size_t eg_size = flux->eg_size;
	size_t xm_size = flux->xm_size;
	size_t ym_size = flux->ym_size;
	size_t zm_size = flux->zm_size;
	int ***cchecker = flux->mapper->cchecker;
	if(g >= eg_size || i >= xm_size || j >= ym_size || k >= zm_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	if(cchecker[k][j][i] & 0b00000001){
		fprintf(stderr, "No mtrl filled in.\n");
		exit(-1);
	}
	#endif
	size_t rt_size = flux->rt_size;
	size_t idx = g*rt_size + mapper_get1Didx(flux->mapper, i, j, k);
	return flux->data[idx];
}

double flux_sumup(const FLUX *flux)
{
	double s = 0.0;
	size_t eg_size = flux->eg_size;
	size_t rt_size = flux->rt_size;
	for(size_t r=0; r<eg_size*rt_size; ++r)
		s += flux->data[r];
	return s;
}
