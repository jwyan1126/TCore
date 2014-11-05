#include"tnsol.h"
#include<stdlib.h>

TNSOL *tnsol_create(size_t eg_size)
{
	TNSOL *tn = malloc(sizeof(TNSOL));
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
	tn->adfgi = calloc(eg_size, sizeof(double));
	tn->adfgj = calloc(eg_size, sizeof(double));
	tn->J = calloc(eg_size, sizeof(double));
	tn->agi1 = calloc(eg_size, sizeof(double));
	tn->agj1 = calloc(eg_size, sizeof(double));
	tn->agi2 = calloc(eg_size, sizeof(double));
	tn->agj2 = calloc(eg_size, sizeof(double));
	tn->agi3 = calloc(eg_size, sizeof(double));
	tn->agj3 = calloc(eg_size, sizeof(double));
	tn->agi4 = calloc(eg_size, sizeof(double));
	tn->agj4 = calloc(eg_size, sizeof(double));
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
	free(tn->adfgi);
	free(tn->adfgj);
	free(tn->J);
	free(tn->agi1);
	free(tn->agj1);
	free(tn->agi2);
	free(tn->agj2);
	free(tn->agi3);
	free(tn->agj3);
	free(tn->agi4);
	free(tn->agj4);
	free(tn);
}

void tnsol_coef_fprintf(const TNSOL *tn, FILE *stream)
{
	size_t eg_size = tn->eg_size;
	for(size_t g=0; g<eg_size; ++g){
		fprintf(stream, "g=%zd\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
			g,tn->agi1[g],tn->agi2[g],tn->agi3[g],tn->agi4[g],tn->agj1[g],tn->agj2[g],tn->agj3[g],tn->agj4[g]);
	}
}
