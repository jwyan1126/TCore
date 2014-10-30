#include"ssol.h"

#include<stdlib.h>

SSOL *ssol_create(MAPPER *mapper)
{
	SSOL *ssol = malloc(sizeof(SSOL));
	size_t eg_size = mapper->eg_size;
	size_t rt_size = mapper->rt_size;
	ssol->eg_size = eg_size;
	ssol->rt_size = rt_size;
	ssol->xm_size = mapper->xm_size;
	ssol->ym_size = mapper->ym_size;
	ssol->zm_size = mapper->zm_size;
	ssol->keff = 1.0;
	ssol->flux = flux_create(mapper);
	for(size_t i=0; i<eg_size * rt_size; ++i)
		ssol->flux->data[i] = 1.0;
	return ssol;
}

void ssol_free(SSOL *ssol)
{
	flux_free(ssol->flux);
	free(ssol);
}

void ssol_fprintf(const SSOL *ssol, FILE *stream)
{
	fprintf(stream, "Keff: %g\n\n", ssol->keff);
	flux_fprintf(ssol->flux, stream);
}
