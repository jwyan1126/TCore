#include"sconf.h"
#include<stddef.h>

int has_xlspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan);
int has_xrspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan);
int has_ylspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan);
int has_yrspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan);
int has_zlspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan);
int has_zrspan(int ***mtrl_set, size_t xspan_size, size_t yspan_size, size_t zspan_size,
		size_t xspan, size_t yspan, size_t zspan);
void bdy_check(int ***bdy_checker, SCONF *sconf);
