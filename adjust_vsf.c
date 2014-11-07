#include"steady_solver.h"

void adjust_vsf(CDAT4 *vsf, double keff)
{
	size_t xm_size = vsf->xsize;
	size_t ym_size = vsf->ysize;
	size_t zm_size = vsf->zsize;
	size_t eg_size = vsf->gsize;
	for(size_t k=0; k<zm_size; ++k)
		for(size_t j=0; j<ym_size; ++j)
			for(size_t i=0; i<xm_size; ++i)
				for(size_t g=0; g<eg_size; ++g)
					vsf->data[k][j][i][g] /= keff;
}
