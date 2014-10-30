#ifndef JCUR_H
#define JCUR_H

typedef struct
{
	size_t eg_size;
	size_t xm_size;
	size_t ym_size;
	size_t zm_size;
	size_t rt_size;
	EDAT4 *jx;
	EDAT4 *jy;
	EDAT4 *jz;
	int ***xchecker;
	int ***ychecker;
	int ***zchecker;
} JCUR;

JCUR *jcur_create(MAPPER *mapper);

void jcur_free(JCUR *jcur);

void cal_jcur(JCUR *jcur, const MESH *mesh, const LEAK *leak, const SSOL *ssol);

#endif
