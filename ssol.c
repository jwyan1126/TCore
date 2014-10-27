#include"solution.h"

SSOL *ssol_create(size_t eg_size, size_t rt_size)
{
	SSOL *ssol = malloc(sizeof(SSOL));
	ssol->keff = 1.0;
	ssol->rt_size = rt_size;
	ssol->eg_size = eg_size;
	ssol->flux = malloc(rt_size * eg_size * sizeof(double));
	for(size_t i=0; i< rt_size * eg_size; ++i)
		ssol->flux[i] = 1.0;
	return ssol;
}

void ssol_free(SSOL *ssol)
{
	free(ssol->flux);
	free(ssol);
}
