#include"checker.h"
#include<stdio.h>
#include<stdlib.h>

// private func.
int has_xlspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan)
{
	if(xspan == 0) return 0;
	if(mtrl_set[xspan-1][yspan][zspan] < 0) return 0;
	return 1;
}

// private func.
int has_xrspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan)
{
	if(xspan == xspan_size-1) return 0;
	if(mtrl_set[xspan+1][yspan][zspan] < 0) return 0;
	return 1;
}

// private func.
int has_ylspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan)
{
	if(yspan == 0) return 0;
	if(mtrl_set[xspan][yspan-1][zspan] < 0) return 0;
	return 1;
}

// private func.
int has_yrspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan)
{
	if(yspan == yspan_size-1) return 0;
	if(mtrl_set[xspan][yspan+1][zspan] < 0) return 0;
	return 1;
}

// private func.
int has_zlspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan)
{
	if(zspan == 0) return 0;
	if(mtrl_set[xspan][yspan][zspan-1] < 0) return 0;
	return 1;
}

// private func.
int has_zrspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan)
{
	if(zspan == zspan_size-1) return 0;
	if(mtrl_set[xspan][yspan][zspan+1] < 0) return 0;
	return 1;
}

// private func.
void bdy_check(int ***bdy_checker, SCONF *sconf)
{
	int ***mtrl_set = sconf->mtrl_set;
	size_t xspan_size = sconf->xm_span_size;
	size_t yspan_size = sconf->ym_span_size;
	size_t zspan_size = sconf->zm_span_size;
	for(size_t zspan=0; zspan<zspan_size; ++zspan)
		for(size_t yspan=0; yspan<yspan_size; ++yspan)
			for(size_t xspan=0; xspan<xspan_size; ++xspan){
				MBLOCK mblock = sconf_get_mblock(sconf, xspan, yspan, zspan);
				for(size_t k = mblock.start_z; k <= mblock.end_z; ++k)
					for(size_t j = mblock.start_y; j <= mblock.end_y; ++j)
						for(size_t i = mblock.start_x; i <= mblock.end_x; ++i){
						int mtrl_id = mtrl_set[xspan][yspan][zspan];
						if(mtrl_id < 0){
							bdy_checker[k][j][i] = 0b00000001;
							continue;
						}
						else
							bdy_checker[k][j][i] = 0b00000010;

						if(i == mblock.start_x && !has_xlspan(mtrl_set,xspan_size,yspan_size,zspan_size,xspan,yspan,zspan))
							bdy_checker[k][j][i] |= 0b00000100;
						if(i == mblock.end_x && !has_xrspan(mtrl_set,xspan_size,yspan_size,zspan_size,xspan,yspan,zspan))
							bdy_checker[k][j][i] |= 0b00001000;
						if(j == mblock.start_y && !has_ylspan(mtrl_set,xspan_size,yspan_size,zspan_size,xspan,yspan,zspan))
							bdy_checker[k][j][i] |= 0b00010000;
						if(j == mblock.end_y && !has_yrspan(mtrl_set,xspan_size,yspan_size,zspan_size,xspan,yspan,zspan))
							bdy_checker[k][j][i] |= 0b00100000;
						if(k == mblock.start_z && !has_zlspan(mtrl_set,xspan_size,yspan_size,zspan_size,xspan,yspan,zspan))
							bdy_checker[k][j][i] |= 0b01000000;
						if(k == mblock.end_z && !has_zrspan(mtrl_set,xspan_size,yspan_size,zspan_size,xspan,yspan,zspan))
							bdy_checker[k][j][i] |= 0b10000000;
						}
			}
}
