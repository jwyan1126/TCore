#include"leak.h"
#include<stdlib.h>
#include<stddef.h>

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
	cdat4_free(leak->lx0);
	cdat4_free(leak->lx1);
	cdat4_free(leak->lx2);
	cdat4_free(leak->ly0);
	cdat4_free(leak->ly1);
	cdat4_free(leak->ly2);
	cdat4_free(leak->lz0);
	cdat4_free(leak->lz1);
	cdat4_free(leak->lz2);
	free(leak);
}

void cal_leakage(LEAK *leak, const MESH *mesh, const EDAT4 *jn)
{
	int xl_bdy = mesh->xl_bdy;
	int xr_bdy = mesh->xr_bdy;
	int yl_bdy = mesh->yl_bdy;
	int yr_bdy = mesh->yr_bdy;
	int zl_bdy = mesh->zl_bdy;
	int zr_bdy = mesh->zr_bdy;
	size_t eg_size = mesh->eg_size;
	size_t xm_size = mesh->xm_mesh_size;
	size_t ym_size = mesh->ym_mesh_size;
	size_t zm_size = mesh->zm_mesh_size;
	int ***cchecker = mesh->cchecker;
	for(size_t k=0; k<zm_size; ++k)
		for(size_t j=0; j<ym_size; ++j)
			for(size_t i=0; i<xm_size; ++i){
				if(cchecker[k][j][i] & 0b00000001) continue;
				double dx = cdat3_get_val(mesh->dx, i, j, k);
				double dy = cdat3_get_val(mesh->dy, i, j, k);
				double dz = cdat3_get_val(mesh->dz, i, j, k);
				for(size_t g=0; g<eg_size; ++g){
					double Jgxl = edat4_get_xlval(jn, g, i, j, k);
					double Jgxr = edat4_get_xrval(jn, g, i, j, k);
					double Jgyl = edat4_get_ylval(jn, g, i, j, k);
					double Jgyr = edat4_get_yrval(jn, g, i, j, k);
					double Jgzl = edat4_get_zlval(jn, g, i, j, k);
					double Jgzr = edat4_get_zrval(jn, g, i, j, k);
					double lgx0 = (Jgyr - Jgyl) / dy + (Jgzr - Jgzl) / dz;
					double lgy0 = (Jgzr - Jgzl) / dz + (Jgxr - Jgxl) / dx;
					double lgz0 = (Jgxr - Jgxl) / dx + (Jgyr - Jgyl) / dy;
					cdat4_set_val(leak->lx0, g, i, j, k, lgx0);
					cdat4_set_val(leak->ly0, g, i, j, k, lgy0);
					cdat4_set_val(leak->lz0, g, i, j, k, lgz0);
				}
			}
	// x
	for(size_t k=0; k<zm_size; ++k)
		for(size_t j=0; j<ym_size; ++j)
			for(size_t i=0; i<xm_size; ++i){
				if(cchecker[k][j][i] & 0b00000001) continue;
				for(size_t g=0; g<eg_size; ++g){
					double Ll, Lm, Lr;
					double dl, dm, dr;
					if(!(cchecker[k][j][i] & 0b00000100)){
						Ll = cdat4_get_val(leak->lx0, g, i-1, j, k);
						dl = cdat3_get_val(mesh->dx, i-1, j, k);
					}
					else{
						if(xl_bdy == 0){
							Ll = cdat4_get_val(leak->lx0, g, i, j, k);
							dl = cdat3_get_val(mesh->dx, i, j, k);
						}
						else{
							Ll = 0.0;
							dl = 0.0;
						}
					}
					if(!(cchecker[k][j][i] & 0b00001000)){
						Lr = cdat4_get_val(leak->lx0, g, i+1, j, k);
						dr = cdat3_get_val(mesh->dx, i+1, j, k);
					}
					else{
						if(xr_bdy == 0){
							Lr = cdat4_get_val(leak->lx0, g, i, j, k);
							dr = cdat3_get_val(mesh->dx, i, j, k);
						}
						else{
							Lr = 0.0;
							dr = 0.0;
						}
					}
					Lm = cdat4_get_val(leak->lx0, g, i, j, k);
					dm = cdat3_get_val(mesh->dx, i, j, k);
					double d = 1.0 / (2.0*(dm+dl)*(dm+dr)*(dm+dl+dr));
					double rou1 = d*dm*((dm+dr)*(dm+2.0*dr)*(Lm-Ll) + (dm+dl)*(dm+2.0*dl)*(Lr-Lm));
					cdat4_set_val(leak->lx1, g, i, j, k, rou1);
					double rou2 = d*dm*dm*((dm+dl)*(Lr-Lm) + (dm+dr)*(Ll-Lm));
					cdat4_set_val(leak->lx2, g, i, j, k, rou2);
				}
			}
	// y
	for(size_t k=0; k<zm_size; ++k)
		for(size_t j=0; j<ym_size; ++j)
			for(size_t i=0; i<xm_size; ++i){
				if(cchecker[k][j][i] & 0b00000001) continue;
				for(size_t g=0; g<eg_size; ++g){
					double Ll, Lm, Lr;
					double dl, dm, dr;
					if(!(cchecker[k][j][i] & 0b00010000)){
						Ll = cdat4_get_val(leak->ly0, g, i, j-1, k);
						dl = cdat3_get_val(mesh->dy, i, j-1, k);
					}
					else{
						if(yl_bdy == 0){
							Ll = cdat4_get_val(leak->ly0, g, i, j, k);
							dl = cdat3_get_val(mesh->dy, i, j, k);
						}
						else{
							Ll = 0.0;
							dl = 0.0;
						}
					}
					if(!(cchecker[k][j][i] & 0b00100000)){
						Lr = cdat4_get_val(leak->ly0, g, i, j+1, k);
						dr = cdat3_get_val(mesh->dy, i, j+1, k);
					}
					else{
						if(yr_bdy == 0){
							Lr = cdat4_get_val(leak->ly0, g, i, j, k);
							dr = cdat3_get_val(mesh->dy, i, j, k);
						}
						else{
							Lr = 0.0;
							dr = 0.0;
						}
					}
					Lm = cdat4_get_val(leak->ly0, g, i, j, k);
					dm = cdat3_get_val(mesh->dy, i, j, k);
					double d = 1.0 / (2.0*(dm+dl)*(dm+dr)*(dm+dl+dr));
					double rou1 = d*dm*((dm+dr)*(dm+2.0*dr)*(Lm-Ll) + (dm+dl)*(dm+2.0*dl)*(Lr-Lm));
					cdat4_set_val(leak->ly1, g, i, j, k, rou1);
					double rou2 = d*dm*dm*((dm+dl)*(Lr-Lm) + (dm+dr)*(Ll-Lm));
					cdat4_set_val(leak->ly2, g, i, j, k, rou2);
				}
			}
	// z
	for(size_t k=0; k<zm_size; ++k)
		for(size_t j=0; j<ym_size; ++j)
			for(size_t i=0; i<xm_size; ++i){
				if(cchecker[k][j][i] & 0b00000001) continue;
				for(size_t g=0; g<eg_size; ++g){
					double Ll, Lm, Lr;
					double dl, dm, dr;
					if(!(cchecker[k][j][i] & 0b01000000)){
						Ll = cdat4_get_val(leak->lz0, g, i, j, k-1);
						dl = cdat3_get_val(mesh->dz, i, j, k-1);
					}
					else{
						if(zl_bdy == 0){
							Ll = cdat4_get_val(leak->lz0, g, i, j, k);
							dl = cdat3_get_val(mesh->dz, i, j, k);
						}
						else{
							Ll = 0.0;
							dl = 0.0;
						}
					}
					if(!(cchecker[k][j][i] & 0b10000000)){
						Lr = cdat4_get_val(leak->lz0, g, i, j, k+1);
						dr = cdat3_get_val(mesh->dz, i, j, k+1);
					}
					else{
						if(zr_bdy == 0){
							Lr = cdat4_get_val(leak->lz0, g, i, j, k);
							dr = cdat3_get_val(mesh->dz, i, j, k);
						}
						else{
							Lr = 0.0;
							dr = 0.0;
						}
					}
					Lm = cdat4_get_val(leak->lz0, g, i, j, k);
					dm = cdat3_get_val(mesh->dz, i, j, k);
					double d = 1.0 / (2.0*(dm+dl)*(dm+dr)*(dm+dl+dr));
					double rou1 = d*dm*((dm+dr)*(dm+2.0*dr)*(Lm-Ll) + (dm+dl)*(dm+2.0*dl)*(Lr-Lm));
					cdat4_set_val(leak->lz1, g, i, j, k, rou1);
					double rou2 = d*dm*dm*((dm+dl)*(Lr-Lm) + (dm+dr)*(Ll-Lm));
					cdat4_set_val(leak->lz2, g, i, j, k, rou2);
				}
			}
}
