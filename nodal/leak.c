#include"leak.h"

LEAK *leak_create(MAPPER *mapper)
{
	LEAK *leak = malloc(sizeof(LEAK));
	leak->eg_size = mapper->eg_size;
	leak->xm_size = mapper->xm_size;
	leak->ym_size = mapper->ym_size;
	leak->zm_size = mapper->zm_size;
	leak->rt_size = mapper->rt_size;
	leak->lx0 = cdat4_create(mapper);
	leak->lx1 = cdat4_create(mapper);
	leak->lx2 = cdat4_create(mapper);
	leak->ly0 = cdat4_create(mapper);
	leak->ly1 = cdat4_create(mapper);
	leak->ly2 = cdat4_create(mapper);
	leak->lz0 = cdat4_create(mapper);
	leak->lz1 = cdat4_create(mapper);
	leak->lz2 = cdat4_create(mapper);
	leak->cchecker = mapper->cchecker;
	return leak;
}

void leak_free(LEAK *leak)
{
	cdat4_free(lx0);
	cdat4_free(lx1);
	cdat4_free(lx2);
	cdat4_free(ly0);
	cdat4_free(ly1);
	cdat4_free(ly2);
	cdat4_free(lz0);
	cdat4_free(lz1);
	cdat4_free(lz2);
	free(leak);
}
