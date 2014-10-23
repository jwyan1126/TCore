#ifndef RECT_MAPPER_H
#define RECT_MAPPER_H

extern size_t EG_SIZE;
extern size_t XM_SIZE;
extern size_t YM_SIZE;
extern size_t ZM_SIZE;
extern size_t RT_SIZE;

typedef struct
{
	size_t xi;
	size_t yi;
	size_t zi;
} XYZ_IDX;

typedef struct
{
	XYZ_IDX *one2three;
	size_t ***three2one; //indexed [k][j][i]
} RECT_MAPPER;

RECT_MAPPER *rect_mapper_create();

void rect_mapper_free(RECT_MAPPER *rect_mapper);

XYZ_IDX rect_mapper_get3Didx(RECT_MAPPER *mapper, size_t idx1D);

size_t rect_mapper_get1Didx(RECT_MAPPER *mapper,size_t i, size_t j, size_t k);


#endif
