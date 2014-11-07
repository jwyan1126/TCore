#ifndef TSOL_H
#define TSOL_H

#include<stddef.h>
#include"flux.h"
#include"pcs.h"
#include"pre_proc/tconf.h"
#include"ssol.h"
#include"pre_proc/mesh.h"

typedef struct
{
	size_t eg_size;
	size_t rt_size;
	size_t xm_size;
	size_t ym_size;
	size_t zm_size;
	size_t pcs_size;
	double time;
	double keff;
	FLUX *flux;
	PCS *pcs;
} TSOL;

TSOL *tsol_create(const TCONF *tconf, const MESH *mesh, const SSOL *ssol);

void tsol_free(TSOL *tsol);

#endif
