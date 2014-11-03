#include"jcur.h"
#include"tnsol.h"
#include"SANM/sanm.h"

void cal_jcur(EDAT4 *jcur, const MESH *mesh, const LEAK *leak, const SSOL *ssol)
{
	size_t eg_size = leak->eg_size;
	size_t xm_size = leak->xm_size;
	size_t ym_size = leak->ym_size;
	size_t zm_size = leak->zm_size;
	int ***xchecker = jcur->xchecker;
	int ***ychecker = jcur->ychecker;
	int ***zchecker = jcur->zchecker;
	TNSOL *tn = tnsol_create(eg_size);
	tn->keff = ssol->keff;
	// x direction [k][j][i]
	for(size_t k=0; k< zm_size; ++k){
		for(size_t j=0; j< ym_size; ++j){
			for(size_t i=0; i< xm_size+1; ++i){
				if(!(xchecker[k][j][i] & 0b0001)) continue;
				if(!(xchecker[k][j][i] & 0b0010)){
					tn->dui = cdat3_get_val(mesh->dx, i-1, j, k);
					for(size_t g=0; g<eg_size; ++g){
						tn->Dgi[g] = cdat4_get_val(mesh->dcoef, g, i-1, j, k);
						tn->vsfgi[g] = cdat4_get_val(mesh->vsf, g, i-1, j, k);
						tn->srgi[g] = cdat4_get_val(mesh->sr, g, i-1, j, k);
						for(size_t from_g=0; from_g<eg_size; ++from_g)
							tn->ssgi[g][from_g] = cdat5_get_val(mesh->ss, g, from_g, i-1, j, k);
						tn->chigi[g] = cdat4_get_val(mesh->chi, g, i-1, j, k);
						tn->phigi[g] = flux_get_val(ssol->flux, g, i-1, j, k);
						tn->lgi0[g] = cdat4_get_val(leak->lx0, g, i-1, j, k);
						tn->lgi1[g] = cdat4_get_val(leak->lx1, g, i-1, j, k);
						tn->lgi2[g] = cdat4_get_val(leak->lx2, g, i-1, j, k);
						tn->adfgi[g] = cdat4_get_val(mesh->adfxr, g, i-1, j, k);
					}
				}
				if(!(xchecker[k][j][i] & 0b0100)){
					tn->duj = cdat3_get_val(mesh->dx, i, j, k);
					for(size_t g=0; g<eg_size; ++g){
						tn->Dgj[g] = cdat4_get_val(mesh->dcoef, g, i, j, k);
						tn->vsfgj[g] = cdat4_get_val(mesh->vsf, g, i, j, k);
						tn->srgj[g] = cdat4_get_val(mesh->sr, g, i, j, k);
						for(size_t from_g=0; from_g<eg_size; ++from_g)
							tn->ssgj[g][from_g] = cdat5_get_val(mesh->ss, g, from_g, i, j, k);
						tn->chigj[g] = cdat4_get_val(mesh->chi, g, i, j, k);
						tn->phigj[g] = flux_get_val(ssol->flux, g, i, j, k);
						tn->lgj0[g] = cdat4_get_val(leak->lx0, g, i, j, k);
						tn->lgj1[g] = cdat4_get_val(leak->lx1, g, i, j, k);
						tn->lgj2[g] = cdat4_get_val(leak->lx2, g, i, j, k);
						tn->adfgj[g] = cdat4_get_val(mesh->adfxl, g, i, j, k);
					}
				}
				if(xchecker[k][j][i] & 0b0010){ tn->bdy = mesh->xl_bdy; sanm_left(tn); }
				else if(xchecker[k][j][i] & 0b0100) { tn->bdy = mesh->xr_bdy; sanm_right(tn); }
				else	sanm_inner(tn);
				for(size_t g=0; g<eg_size; ++g)
					jcur->xdata[k][j][i][g] = tn->J[g];
			}
		}
	}
	// y direction [i][k][j]
	for(size_t i=0; i< xm_size; ++i){
		for(size_t k=0; k< zm_size; ++k){
			for(size_t j=0; j< ym_size+1; ++j){
				if(!(ychecker[i][k][j] & 0b0001)) continue;
				if(!(ychecker[i][k][j] & 0b0010)){
					tn->dui = cdat3_get_val(mesh->dy, i, j-1, k);
					for(size_t g=0; g<eg_size; ++g){
						tn->Dgi[g] = cdat4_get_val(mesh->dcoef, g, i, j-1, k);
						tn->vsfgi[g] = cdat4_get_val(mesh->vsf, g, i, j-1, k);
						tn->srgi[g] = cdat4_get_val(mesh->sr, g, i, j-1, k);
						for(size_t from_g=0; from_g<eg_size; ++from_g)
							tn->ssgi[g][from_g] = cdat5_get_val(mesh->ss, g, from_g, i, j-1, k);
						tn->chigi[g] = cdat4_get_val(mesh->chi, g, i, j-1, k);
						tn->phigi[g] = flux_get_val(ssol->flux, g, i, j-1, k);
						tn->lgi0[g] = cdat4_get_val(leak->lx0, g, i, j-1, k);
						tn->lgi1[g] = cdat4_get_val(leak->lx1, g, i, j-1, k);
						tn->lgi2[g] = cdat4_get_val(leak->lx2, g, i, j-1, k);
						tn->adfgi[g] = cdat4_get_val(mesh->adfyr, g, i, j-1, k);
					}
				}
				if(!(ychecker[i][k][j] & 0b0100)){
					tn->duj = cdat3_get_val(mesh->dy, i, j, k);
					for(size_t g=0; g<eg_size; ++g){
						tn->Dgj[g] = cdat4_get_val(mesh->dcoef, g, i, j, k);
						tn->vsfgj[g] = cdat4_get_val(mesh->vsf, g, i, j, k);
						tn->srgj[g] = cdat4_get_val(mesh->sr, g, i, j, k);
						for(size_t from_g=0; from_g<eg_size; ++from_g)
							tn->ssgj[g][from_g] = cdat5_get_val(mesh->ss, g, from_g, i, j, k);
						tn->chigj[g] = cdat4_get_val(mesh->chi, g, i, j, k);
						tn->phigj[g] = flux_get_val(ssol->flux, g, i, j, k);
						tn->lgj0[g] = cdat4_get_val(leak->lx0, g, i, j, k);
						tn->lgj1[g] = cdat4_get_val(leak->lx1, g, i, j, k);
						tn->lgj2[g] = cdat4_get_val(leak->lx2, g, i, j, k);
						tn->adfgj[g] = cdat4_get_val(mesh->adfyl, g, i, j, k);
					}
				}
				if(ychecker[i][k][j] & 0b0010) { tn->bdy = mesh->yl_bdy; sanm_left(tn); }
				else if(ychecker[i][k][j] & 0b0100) { tn->bdy = mesh->yr_bdy; sanm_right(tn); }
				else sanm_inner(tn);
				for(size_t g=0; g<eg_size; ++g)
					jcur->ydata[i][k][j][g] = tn->J[g];
			}
		}
	}
	// z direction [j][i][k]
	for(size_t j=0; j< ym_size; ++j){
		for(size_t i=0; i< xm_size; ++i){
			for(size_t k=0; k< zm_size+1; ++k){
				if(!(zchecker[j][i][k] & 0b0001)) continue;
				if(!(zchecker[j][i][k] & 0b0010)){
					tn->dui = cdat3_get_val(mesh->dz, i, j, k-1);
					for(size_t g=0; g<eg_size; ++g){
						tn->Dgi[g] = cdat4_get_val(mesh->dcoef, g, i, j, k-1);
						tn->vsfgi[g] = cdat4_get_val(mesh->vsf, g, i, j, k-1);
						tn->srgi[g] = cdat4_get_val(mesh->sr, g, i, j, k-1);
						for(size_t from_g=0; from_g<eg_size; ++from_g)
							tn->ssgi[g][from_g] = cdat5_get_val(mesh->ss, g, from_g, i, j, k-1);
						tn->chigi[g] = cdat4_get_val(mesh->chi, g, i, j, k-1);
						tn->phigi[g] = flux_get_val(ssol->flux, g, i, j, k-1);
						tn->lgi0[g] = cdat4_get_val(leak->lx0, g, i, j, k-1);
						tn->lgi1[g] = cdat4_get_val(leak->lx1, g, i, j, k-1);
						tn->lgi2[g] = cdat4_get_val(leak->lx2, g, i, j, k-1);
						tn->adfgi[g] = cdat4_get_val(mesh->adfzr, g, i, j, k-1);
					}
				}
				if(!(zchecker[j][i][k] & 0b0100)){
					tn->duj = cdat3_get_val(mesh->dz, i, j, k);
					for(size_t g=0; g<eg_size; ++g){
						tn->Dgj[g] = cdat4_get_val(mesh->dcoef, g, i, j, k);
						tn->vsfgj[g] = cdat4_get_val(mesh->vsf, g, i, j, k);
						tn->srgj[g] = cdat4_get_val(mesh->sr, g, i, j, k);
						for(size_t from_g=0; from_g<eg_size; ++from_g)
							tn->ssgj[g][from_g] = cdat5_get_val(mesh->ss, g, from_g, i, j, k);
						tn->chigj[g] = cdat4_get_val(mesh->chi, g, i, j, k);
						tn->phigj[g] = flux_get_val(ssol->flux, g, i, j, k);
						tn->lgj0[g] = cdat4_get_val(leak->lx0, g, i, j, k);
						tn->lgj1[g] = cdat4_get_val(leak->lx1, g, i, j, k);
						tn->lgj2[g] = cdat4_get_val(leak->lx2, g, i, j, k);
						tn->adfgj[g] = cdat4_get_val(mesh->adfzl, g, i, j, k);
					}
				}
				if(zchecker[j][i][k] & 0b0010) { tn->bdy = mesh->zl_bdy; sanm_left(tn); }
				else if(zchecker[j][i][k] & 0b0100) { tn->bdy = mesh->zr_bdy; sanm_right(tn); }
				else sanm_inner(tn);
				for(size_t g=0; g<eg_size; ++g)
					jcur->zdata[j][i][k][g] = tn->J[g];
			}
		}
	}
	tnsol_free(tn);
}
