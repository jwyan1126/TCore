#include"steady_solver.h"
#include<stdlib.h>

void cal_DFDM(EDAT4 *DFDM, const SCONF *sconf, const MESH *mesh)
{
	size_t eg_size = sconf->eg_size;
	size_t xm_size = sconf->xm_mesh_size;
	size_t ym_size = sconf->ym_mesh_size;
	size_t zm_size = sconf->zm_mesh_size;
	int xl_bdy = sconf->xl_bdy;
	int xr_bdy = sconf->xr_bdy;
	int yl_bdy = sconf->yl_bdy;
	int yr_bdy = sconf->yr_bdy;
	int zl_bdy = sconf->zl_bdy;
	int zr_bdy = sconf->zr_bdy;
	CDAT3 *dx = mesh->dx;
	CDAT3 *dy = mesh->dy;
	CDAT3 *dz = mesh->dz;
	CDAT4 *dcoef = mesh->dcoef;
	int ***xchecker = DFDM->xchecker;
	int ***ychecker = DFDM->ychecker;
	int ***zchecker = DFDM->zchecker;
	for(size_t k=0; k<zm_size; ++k)
		for(size_t j=0; j<ym_size; ++j)
			for(size_t i=0; i<xm_size+1; ++i){
				//xdata[k][j][i]
				if(!(xchecker[k][j][i] & 0b0001)) continue;
				for(size_t g=0; g<eg_size; ++g){
					if(xchecker[k][j][i] & 0b0010){
						if(xchecker[k][j][i] & 0b0100) {fprintf(stderr, "inner error\n"); exit(-1);}//DEBUG
						if(xl_bdy == 0)
							DFDM->xdata[k][j][i][g] = 0.0;
						else if(xl_bdy == 1){
							double beta_r = cdat4_get_val(dcoef, g, i, j, k) / cdat3_get_val(dx, i, j, k);
							DFDM->xdata[k][j][i][g] = 2.0 * beta_r;
						}
						else if(xl_bdy == 2){
							double beta_r = cdat4_get_val(dcoef, g, i, j, k) / cdat3_get_val(dx, i, j, k);
							DFDM->xdata[k][j][i][g] = 0.5 * beta_r / (0.25 + beta_r);
						}
						else {fprintf(stderr, "bdy unrecognized.\n"); exit(-1);}
					}
					else if(xchecker[k][j][i] & 0b0100){
						if(xchecker[k][j][i] & 0b0010) {fprintf(stderr, "inner error\n"); exit(-1);}
						if(xr_bdy == 0)
							DFDM->xdata[k][j][i][g] = 0.0;
						else if(xr_bdy == 1){
							double beta_l = cdat4_get_val(dcoef, g, i-1, j, k) / cdat3_get_val(dx, i-1, j, k);
							DFDM->xdata[k][j][i][g] = 2.0 * beta_l;
						}
						else if(xr_bdy == 2){
							double beta_l = cdat4_get_val(dcoef, g, i-1, j, k) / cdat3_get_val(dx, i-1, j, k);
							DFDM->xdata[k][j][i][g] = 0.5 * beta_l / (0.25 + beta_l);
						}
						else {fprintf(stderr, "bdy unrecognized.\n"); exit(-1);}
					}
					else{
						double beta_l = cdat4_get_val(dcoef, g, i-1, j, k) / cdat3_get_val(dx, i-1, j, k);
						double beta_r = cdat4_get_val(dcoef, g, i, j, k) / cdat3_get_val(dx, i, j, k);
						DFDM->xdata[k][j][i][g] = 2.0 * beta_l * beta_r / (beta_l + beta_r);
					}
				}
			}
	for(size_t i=0; i<xm_size; ++i)
		for(size_t k=0; k<zm_size; ++k)
			for(size_t j=0; j<ym_size+1; ++j){
				//ydata[i][k][j]
				if(!(ychecker[i][k][j] & 0b0001)) continue;
				for(size_t g=0; g<eg_size; ++g){
					if(ychecker[i][k][j] & 0b0010){
						if(yl_bdy == 0)
							DFDM->ydata[i][k][j][g] = 0.0;
						else if(yl_bdy == 1){
							double beta_r = cdat4_get_val(dcoef, g, i, j, k) / cdat3_get_val(dy, i, j, k);
							DFDM->ydata[i][k][j][g] = 2.0 * beta_r;
						}
						else if(yl_bdy == 2){
							double beta_r = cdat4_get_val(dcoef, g, i, j, k) / cdat3_get_val(dy, i, j, k);
							DFDM->ydata[i][k][j][g] = 0.5 * beta_r / (0.25 + beta_r);
						}
						else {fprintf(stderr, "bdy unrecognized.\n"); exit(-1);}
					}
					else if(ychecker[i][k][j] & 0b0100){
						if(yr_bdy == 0)
							DFDM->ydata[i][k][j][g] = 0.0;
						else if(yr_bdy == 1){
							double beta_l = cdat4_get_val(dcoef, g, i, j-1, k) / cdat3_get_val(dy, i, j-1, k);
							DFDM->ydata[i][k][j][g] = 2.0 * beta_l;
						}
						else if(yr_bdy == 2){
							double beta_l = cdat4_get_val(dcoef, g, i, j-1, k) / cdat3_get_val(dy, i, j-1, k);
							DFDM->ydata[i][k][j][g] = 0.5 * beta_l / (0.25 + beta_l);
						}
						else {fprintf(stderr, "bdy unrecognized.\n"); exit(-1);}
					}
					else{
						double beta_l = cdat4_get_val(dcoef, g, i, j-1, k) / cdat3_get_val(dy, i, j-1, k);
						double beta_r = cdat4_get_val(dcoef, g, i, j, k) / cdat3_get_val(dy, i, j, k);
						DFDM->ydata[i][k][j][g] = 2.0*beta_l*beta_r/(beta_l + beta_r);
					}
				}
			}
	for(size_t j=0; j<ym_size; ++j)
		for(size_t i=0; i<xm_size; ++i)
			for(size_t k=0; k<zm_size+1; ++k){
				//zdata[j][i][k]
				if(!(zchecker[j][i][k] & 0b0001)) continue;
				for(size_t g=0; g<eg_size; ++g){
					if(zchecker[j][i][k] & 0b0010){
						if(zl_bdy == 0)
							DFDM->zdata[j][i][k][g] = 0.0;
						else if(zl_bdy == 1){
							double beta_r = cdat4_get_val(dcoef, g, i, j, k) / cdat3_get_val(dz, i, j, k);
							DFDM->zdata[j][i][k][g] = 2.0 * beta_r;
						}
						else if(zl_bdy == 2){
							double beta_r = cdat4_get_val(dcoef, g, i, j, k) / cdat3_get_val(dz, i, j, k);
							DFDM->zdata[j][i][k][g] = 0.5 * beta_r / (0.25 + beta_r);
						}
						else {fprintf(stderr, "bdy unrecognized.\n"); exit(-1);}
					}
					else if(zchecker[j][i][k] & 0b0100){
						if(zr_bdy == 0)
							DFDM->zdata[j][i][k][g] = 0.0;
						else if(zr_bdy == 1){
							double beta_l = cdat4_get_val(dcoef, g, i, j, k-1) / cdat3_get_val(dz, i, j, k-1);
							DFDM->zdata[j][i][k][g] = 2.0 * beta_l;
						}
						else if(zr_bdy == 2){
							double beta_l = cdat4_get_val(dcoef, g, i, j, k-1) / cdat3_get_val(dz, i, j, k-1);
							DFDM->zdata[j][i][k][g] = 0.5 * beta_l / (0.25 + beta_l);
						}
						else {fprintf(stderr, "bdy unrecognized.\n"); exit(-1);}
					}
					else{
						double beta_l = cdat4_get_val(dcoef, g, i, j, k-1) / cdat3_get_val(dz, i, j, k-1);
						double beta_r = cdat4_get_val(dcoef, g, i, j, k) / cdat3_get_val(dz, i, j, k);
						DFDM->zdata[j][i][k][g] = 2.0*beta_l*beta_r/(beta_l + beta_r);
					}
				}
			}
}
