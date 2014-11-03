#include"../tnsol.h"
#include"cal_bgk.h"
#include"sanm.h"

int main()
{
	TNSOL *tn = tnsol_create(2);
	tn->eg_size = 2;
	tn->keff = 0.77;
	tn->Dgi[0] = 1.5450;
	tn->Dgi[1] = 0.3126;
	tn->Dgj[0] = 1.5450;
	tn->Dgj[1] = 0.3126;
	tn->dui = 4.0;
	tn->duj = 4.0;
	tn->vsfgi[0] = 0.04;
	tn->vsfgi[1] = 0.1;
	tn->vsfgj[0] = 0.05;
	tn->vsfgj[1] = 0.15;
	tn->phigi[0] = 0.01838;
	tn->phigi[1] = 0.02127;
	tn->phigj[0] = 0.0149;
	tn->phigj[1] = 0.0172;
	tn->srgi[0] = 0.028824;
	tn->srgi[1] = 0.008736;
	tn->srgj[0] = 0.028824;
	tn->srgj[1] = 0.008736;
	tn->ssgi[1][0] = 0.02838;
	tn->ssgj[1][0] = 0.02838;
	tn->chigi[0] = 1.0;
	tn->chigi[1] = 0.0;
	tn->chigj[0] = 1.0;
	tn->chigj[1] = 0.0;

	tn->lgi0[0] = 0.01;
	tn->lgi0[1] = 0.015;
	tn->lgj0[0] = 0.01;
	tn->lgj0[1] = 0.015;

	tn->lgi1[0] = 0.001;
	tn->lgi1[1] = 0.002;
	tn->lgj1[0] = 0.001;
	tn->lgj1[1] = 0.002;

	tn->lgi2[0] = 0.001;
	tn->lgi2[1] = 0.001;
	tn->lgj2[0] = 0.001;
	tn->lgj2[1] = 0.001;

	tn->adfgi[0] = 1.0;
	tn->adfgi[1] = 1.0;
	tn->adfgj[0] = 1.0;
	tn->adfgj[1] = 1.0;

	tn->bdy = 1;

	sanm_right(tn);
	tnsol_free(tn);
}
