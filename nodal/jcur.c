JCUR *jcur_create(MAPPER *mapper)
{
	JCUR *jcur = malloc(sizeof(JCUR));
	jcur->eg_size = mapper->eg_size;
	jcur->xm_size = mapper->xm_size;
	jcur->ym_size = mapper->ym_size;
	jcur->zm_size = mapper->zm_size;
	jcur->rt_size = mapper->rt_size;
	jcur->jx = edat4_create(mapper);
	jcur->jy = edat4_create(mapper);
	jcur->jz = edat4_create(mapper);
	return jcur;
}

void jcur_free(JCUR *jcur)
{
	edat4_free(jcur->jx);
	edat4_free(jcur->jy);
	edat4_free(jcur->jz);
	free(jcur);
}

void cal_jcur(JCUR *jcur, const MESH *mesh, const LEAK *leak, const SSOL *ssol)
{
	size_t eg_size = leak->eg_size;
	size_t xm_size = leak->xm_size;
	size_t ym_size = leak->ym_size;
	size_t zm_size = leak->zm_size;
	int ***xchecker = jcur->xchecker;
	int ***ychecker = jcur->ychecker;
	int ***zchecker = jcur->zchecker;
	TNINP *tn = tninp_create(eg_size);
	// x direction [k][j][i]
	for(size_t k=0; k< zm_size; ++k){
		for(size_t j=0; j< ym_size; ++j){
			for(size_t i=0; i< xm_size+1; ++i){
				if(!xchecker[k][j][i] & 0b0001) continue;
				if(!xchecker[k][j][i] & 0b0010){
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
						tn->lgi1[g] = cdat4_get_val(leak->lx1 g, i-1, j, k);
						tn->lgi2[g] = cdat4_get_val(leak->lx2, g, i-1, j, k);
					}
				}
				if(!xchecker[k][j][i] & 0b0100){
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
						tn->lgj1[g] = cdat4_get_val(leak->lx1 g, i, j, k);
						tn->lgj2[g] = cdat4_get_val(leak->lx2, g, i, j, k);
					}
				}
			}
		}
	}
	// y direction [i][k][j]
	for(size_t i=0; i< xm_size; ++i){
		for(size_t k=0; k< zm_size; ++k){
			for(size_t j=0; j< ym_size+1; ++j){
				if(!ychecker[i][k][j] & 0b0001) continue;
				if(!ychecker[i][k][j] & 0b0010){
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
						tn->lgi1[g] = cdat4_get_val(leak->lx1 g, i, j-1, k);
						tn->lgi2[g] = cdat4_get_val(leak->lx2, g, i, j-1, k);
					}
				}
				if(!ychecker[i][k][j] & 0b0100){

				}
			}
		}
	}
	// z direction [j][i][k]
	for(size_t j=0; j< ym_size; ++j){
		for(size_t i=0; i< xm_size; ++i){
			for(size_t k=0; k< zm_size+1; ++k){
				if(!zchecker[j][i][k] & 0b0001) continue;
				if(!zchecker[j][i][k] & 0b0010){

				}
				if(!zchecker[j][i][k] & 0b0100){

				}
			}
		}
	}
	
	
	tninp_free(tn);
}
