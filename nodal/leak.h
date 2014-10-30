#ifndef LEAK_H
#define LEAK_H

typedef struct
{
	size_t eg_size;
	size_t xm_size;
	size_t ym_size;
	size_t zm_size;
	size_t rt_size;
	CDAT4 *lx0;
	CDAT4 *lx1;
	CDAT4 *lx2;
	CDAT4 *ly0;
	CDAT4 *ly1;
	CDAT4 *ly2;
	CDAT4 *lz0;
	CDAT4 *lz1;
	CDAT4 *lz2;
	int ***cchecker;
} LEAK;

LEAK *leak_create(MAPPER *mapper);

void leak_free(LEAK *leak);

#endif
