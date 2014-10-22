#include<stdio.h>
#include<stdlib.h>
#include"mtrl.h"
#include"mtrllib.h"

size_t EGSIZE;

int main()
{
	EGSIZE = 2;
	double *chi = calloc(EGSIZE, sizeof(double));
	double *dcoef = calloc(EGSIZE, sizeof(double));
	double *sa = calloc(EGSIZE, sizeof(double));
	double *vsf = calloc(EGSIZE, sizeof(double));
	double **ss = malloc(EGSIZE * sizeof(double *));
	for(size_t g=0; g<EGSIZE; ++g)
		ss[g] = calloc(EGSIZE, sizeof(double));
	MTRL *mtrl1 = mtrl_create(1,chi,dcoef,sa,vsf,ss);
	MTRL *mtrl2 = mtrl_create(2,chi,dcoef,sa,vsf,ss);
	MTRLLIB *mlib = mtrllib_create();
	printf("Before Add, mtrllib_size = %zd\n", mtrllib_get_size(mlib));
	mtrllib_add(mlib, mtrl1);
	mtrllib_add(mlib, mtrl2);
	mtrllib_remove_fromid(mlib, 3);
	MTRL *m = mtrllib_get_fromid(mlib, 1);
	if(m)
		printf("Has element\n");
	else
		printf("Hasn't element\n");
	printf("After Add, mtrllib_size = %zd\n", mtrllib_get_size(mlib));
	mtrllib_free(mlib);
	mtrl_free(mtrl1);
	mtrl_free(mtrl2);
	return 0;
}
