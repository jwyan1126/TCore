#ifndef PCS_H
#define PCS_H

#include<stddef.h>
#include<stdio.h>
#include"pre_proc/mapper.h"
#include"pre_proc/mesh.h"
#include"pre_proc/tconf.h"
#include"flux.h"

typedef struct
{
	size_t pcs_size;
	size_t xm_size;
	size_t ym_size;
	size_t zm_size;
	size_t rt_size;
	double **data;
	MAPPER *mapper;
} PCS;

PCS *pcs_create(size_t pcs_size, MAPPER *mapper);

void pcs_free(PCS *pcs);

double pcs_get_val(const PCS *pcs, size_t p, size_t i, size_t j, size_t k);

void pcs_copy(PCS *tar_pcs, const PCS *src_pcs);

void pcs_init(PCS *pcs, FLUX *flux, MESH *mesh, TCONF *tconf);

void pcs_fprintf(const PCS *pcs, FILE *stream);

#endif
