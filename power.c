#include"power.h"
#include<stdlib.h>

POWER *power_create(MAPPER *mapper)
{
	POWER *power = malloc(sizeof(POWER));
	size_t rt_size = mapper->rt_size;
	power->xm_size = mapper->xm_size;
	power->ym_size = mapper->ym_size;
	power->zm_size = mapper->zm_size;
	power->rt_size = rt_size;
	power->data = calloc(rt_size, sizeof(double));
	power->mapper = mapper;
	return power;
}

void power_free(POWER *power)
{
	free(power->data);
	free(power);
}

void power_fprintf(const POWER *power, FILE *stream)
{
	MAPPER *mapper = power->mapper;
	size_t xm_size = mapper->xm_size;
	size_t ym_size = mapper->ym_size;
	size_t zm_size = mapper->zm_size;
	int ***cchecker = mapper->cchecker;
	for(size_t k=0; k<zm_size; ++k){
		fprintf(stream, "k=%zd\n", k);
		for(size_t j=0; j<ym_size; ++j){
			for(size_t i=0; i<xm_size; ++i){
				size_t idx = mapper_get1Didx(mapper, i, j, k);
				double val = (cchecker[k][j][i] & 0b00000001) ? 0.0 : power->data[idx];
				fprintf(stream, "%g\t", val);
			}
			fprintf(stream, "\n");
		}
		fprintf(stream, "\n");
	}
}

void power_normalize(POWER *power)
{
	size_t rt_size = power->rt_size;
	double s = 0.0;
	for(size_t i=0; i<rt_size; ++i)
		s += power->data[i];
	s /= rt_size;
	for(size_t i=0; i<rt_size; ++i)
		power->data[i] /= s;
}

double power_get_val(const POWER *power, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	size_t xm_size = power->xm_size;
	size_t ym_size = power->ym_size;
	size_t zm_size = power->zm_size;
	int ***cchecker = power->mapper->cchecker;
	if(i >= xm_size || j >= ym_size || k >= zm_size){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	if(cchecker[k][j][i] & 0b00000001){
		fprintf(stderr, "No mtrl filled in.\n");
		exit(-1);
	}
	#endif
	size_t idx = mapper_get1Didx(power->mapper, i, j, k);
	return power->data[idx];
}

void power_cal(POWER *power, const FLUX *flux, const CDAT4 *vsf)
{
	MAPPER *mapper = flux->mapper;
	size_t rt_size = mapper->rt_size;
	size_t eg_size = mapper->eg_size;
	for(size_t r=0; r<rt_size; ++r){
		XYZ_IDX xyz = mapper_get3Didx(mapper,r);
		size_t i=xyz.xi; size_t j=xyz.yi; size_t k=xyz.zi;
		double s = 0.0;
		for(size_t g=0; g<eg_size; ++g)
			s += flux->data[g*rt_size+r] * cdat4_get_val(vsf,g,i,j,k);
		power->data[r] = s;
	}
}

double power_sumup(const POWER *power, const CDAT3 *dx, const CDAT3 *dy, const CDAT3 *dz)
{
	MAPPER *mapper = power->mapper;
	size_t rt_size = mapper->rt_size;
	double s = 0.0;
	for(size_t r=0; r<rt_size; ++r){
		XYZ_IDX xyz = mapper_get3Didx(mapper,r);
		size_t i=xyz.xi; size_t j=xyz.yi; size_t k=xyz.zi;
		s += power->data[r] * cdat3_get_val(dx,i,j,k)
				    * cdat3_get_val(dy,i,j,k)
				    * cdat3_get_val(dz,i,j,k);
	}
	return s;
}
