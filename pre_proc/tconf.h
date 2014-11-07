#ifndef TCONF_H
#define TCONF_H

#include<stddef.h>
#include<stdio.h>
#include"input.h"

typedef struct
{
	size_t eg_size;
	size_t pcs_size;
	double *nvel;
	double *lambdas;
	double *betas;
	double beta;
	double tau;
	int steps;
} TCONF;

TCONF *tconf_create(const INPUT *input);

void tconf_free(TCONF *tconf);

void tconf_fprintf(const TCONF *tconf, FILE *stream);

#endif
