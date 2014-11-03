#include"steady_solver.h"

void cal_DNOD(EDAT4 *DNOD, const MESH *mesh, const EDAT4 *DFDM, const EDAT4 *Jn, const SSOL *ssol)
{
	size_t eg_size = DNOD->gsize;
	size_t xm_size = DNOD->xsize;
	size_t ym_size = DNOD->ysize;
	size_t zm_size = DNOD->zsize;
	int ***xchecker = DNOD->xchecker;
	int ***ychecker = DNOD->ychecker;
	int ***zchecker = DNOD->zchecker;
	FLUX *flux = ssol->flux;
	CDAT4 *adfxl = mesh->adfxl;
	CDAT4 *adfxr = mesh->adfxr;
	CDAT4 *adfyl = mesh->adfyl;
	CDAT4 *adfyr = mesh->adfyr;
	CDAT4 *adfzl = mesh->adfzl;
	CDAT4 *adfzr = mesh->adfzr;
	
	// x direction [k][j][i]
	for(size_t k=0; k< zm_size; ++k){
		for(size_t j=0; j< ym_size; ++j){
			for(size_t i=0; i< xm_size+1; ++i){
				if(!xchecker[k][j][i] & 0b0001) continue;
				for(size_t g=0; g<eg_size; ++g){
					double lval, rval;
					if(xchecker[k][j][i] & 0b0010)	lval = 0.0;
					else	lval = flux_get_val(flux, g, i-1, j, k) * cdat4_get_val(adfxr, g, i-1, j, k);
					if(xchecker[k][j][i] & 0b0100)	rval = 0.0;
					else	rval = flux_get_val(flux, g, i, j, k) * cdat4_get_val(adfxl, g, i, j, k);
					double DFDM_val = DFDM->xdata[k][j][i][g];
					double Jn_val = Jn->xdata[k][j][i][g];
					DNOD->xdata[k][j][i][g] = -(DFDM_val * (rval - lval) + Jn_val) / (rval + lval);
				}
			}
		}
	}
	// y direction [i][k][j]
	for(size_t i=0; i< xm_size; ++i){
		for(size_t k=0; k< zm_size; ++k){
			for(size_t j=0; j< ym_size+1; ++j){
				if(!ychecker[i][k][j] & 0b0001) continue;
				for(size_t g=0; g<eg_size; ++g){
					double lval, rval;
					if(ychecker[i][k][j] & 0b0010) lval = 0.0;
					else	lval = flux_get_val(flux, g, i, j-1, k) * cdat4_get_val(adfyr, g, i, j-1, k);
					if(ychecker[i][k][j] & 0b0100) rval = 0.0;
					else	rval = flux_get_val(flux, g, i, j, k) * cdat4_get_val(adfyl, g, i, j, k);
					double DFDM_val = DFDM->ydata[i][k][j][g];
					double Jn_val = Jn->ydata[i][k][j][g];
					DNOD->ydata[i][k][j][g] = -(DFDM_val * (rval - lval) + Jn_val) / (rval + lval);
				}
			}
		}
	}
	// z direction [j][i][k]
	for(size_t j=0; j< ym_size; ++j){
		for(size_t i=0; i< xm_size; ++i){
			for(size_t k=0; k< zm_size+1; ++k){
				if(!zchecker[j][i][k] & 0b0001) continue;
				for(size_t g=0; g<eg_size; ++g){
					double lval, rval;
					if(zchecker[j][i][k] & 0b0010) lval = 0.0;
					else	lval = flux_get_val(flux, g, i, j, k-1) * cdat4_get_val(adfzr, g, i, j, k-1);
					if(zchecker[j][i][k] & 0b0100) rval = 0.0;
					else	rval = flux_get_val(flux, g, i, j, k) * cdat4_get_val(adfzl, g, i, j, k);
					double DFDM_val = DFDM->zdata[j][i][k][g];
					double Jn_val = Jn->zdata[j][i][k][g];
					DNOD->zdata[j][i][k][g] = -(DFDM_val * (rval - lval) + Jn_val) / (rval + lval);
				}
			}
		}
	}
}
