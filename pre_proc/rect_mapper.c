#include"rect_mapper.h"
#include<stdlib.h>

RECT_MAPPER *rect_mapper_create()
{
	RECT_MAPPER *mapper = malloc(sizeof(RECT_MAPPER));
	mapper->one2three = malloc(RT_SIZE*sizeof(XYZ_IDX));
	mapper->three2one = malloc(ZM_SIZE*sizeof(**size_t))
	for(size_t k=0; k<ZM_SIZE; ++k){
		mapper->one2three[k] = malloc(YM_SIZE*sizeof(*size_t))
		for(size_t j=0; j<YM_SIZE; ++j)
			mapper->one2three[k][j] = calloc(XM_SIZE,sizeof(size_t))
	}
	return mapper;
}

void rect_mapper_free(RECT_MAPPER *mapper)
{
	for(size_t k=0; k<ZM_SIZE; ++k){
		for(size_t j=0; j<YM_SIZE; ++j)
			free(mapper->one2three[k][j]);
		free(mapper->one2three[k]);
	}
	free(mapper->one2three);
	free(mapper);
}

inline XYZ_IDX rect_mapper_get3Didx(RECT_MAPPER *mapper, size_t idx1D)
{
	#ifdef DEBUG
	if(idx1D >= RT_SIZE){
		fprintf(stderr, "'idx1D' out of range.\n");
		exit(-1);
	}
	#endif
	return mapper->one2three[idx1D];
}

inline size_t rect_mapper_get1Didx(RECT_MAPPER *mapper, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	if(i >= XM_SIZE || j >= YM_SIZE || k >= ZM_SIZE){
		fprintf(stderr, "'idx3D' out of range.\n");
		exit(-1);
	}
	#endif
	size_t r = three2one[idx.zi][idx.yi][idx.xi];
	#ifdef DEBUG
	if(r <0 || r >=RT_SIZE){
		fprintf(stderr, "'idx3D' out of range.\n");
		exit(-1);
	}
	#endif
	return r;
}
