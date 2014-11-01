#include"tnsol.h"

TNSOL *tnsol_create(size_t eg_size)
{
	TNSOL *tn = malloc(TNSOL);
	tn->eg_size = eg_size;
	tn->Dgi = calloc(eg_size, sizeof(double));
	tn->Dgj = calloc(eg_size, sizeof(double));
	tn->vsfgi = calloc(eg_size, sizeof(double));
	tn->vsfgj = calloc(eg_size, sizeof(double));
	tn->phigi = calloc(eg_size, sizeof(double));
	tn->phigj = calloc(eg_size, sizeof(double));
	tn->srgi = calloc(eg_size, sizeof(double));
	tn->srgj = calloc(eg_size, sizeof(double));
	tn->ssgi = calloc(eg_size, sizeof(double *));
	tn->ssgj = calloc(eg_size, sizeof(double *));
	for(size_t g=0; g<eg_size; ++g){
		tn->ssgi[g] = calloc(eg_size, sizeof(double));
		tn->ssgj[g] = calloc(eg_size, sizeof(double));
	}
	tn->chigi = calloc(eg_size, sizeof(double));
	tn->chigj = calloc(eg_size, sizeof(double));
	tn->lgi0 = calloc(eg_size, sizeof(double));
	tn->lgj0 = calloc(eg_size, sizeof(double));
	tn->lgi1 = calloc(eg_size, sizeof(double));
	tn->lgj1 = calloc(eg_size, sizeof(double));
	tn->lgi2 = calloc(eg_size, sizeof(double));
	tn->lgj2 = calloc(eg_size, sizeof(double));
	tn->J = calloc(eg_size, sizeof(double));
	return tn;
}

void tnsol_free(TNSOL *tn)
{
	size_t eg_size = tn->eg_size;
	free(tn->Dgi);
	free(tn->Dgj);
	free(tn->vsfgi);
	free(tn->vsfgj);
	free(tn->phigi);
	free(tn->phigj);
	free(tn->srgi);
	free(tn->srgj);
	for(size_t g=0; g<eg_size; ++g){
		free(tn->ssgi[g]);
		free(tn->ssgj[g]);
	}
	free(tn->chigi);
	free(tn->chigj);
	free(tn->lgi0);
	free(tn->lgj0);
	free(tn->lgi1);
	free(tn->lgj1);
	free(tn->lgi2);
	free(tn->lgj2);
	free(tn->J);
	free(tn);
}
