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
	flux->mapper = mapper;
	return flux;
}

void flux_free(FLUX *flux)
{
	free(flux->data);
	free(flux);
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
