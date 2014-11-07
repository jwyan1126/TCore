#include"tconf.h"
#include<stdlib.h>

TCONF *tconf_create(const INPUT *input)
{
	TCONF *tconf = malloc(sizeof(TCONF));
	tconf->eg_size = input->eg_size;
	tconf->pcs_size = input->pcs_size;
	tconf->nvel = input->nvel;
	tconf->lambdas = input->lambdas;
	tconf->betas = input->betas;
	tconf->tau = input->tau;
	tconf->steps = input->steps;
	tconf->beta = 0;
	for(size_t p=0; p<input->pcs_size; ++p)
		tconf->beta += tconf->betas[p];
	return tconf;
}

void tconf_free(TCONF *tconf)
{
	free(tconf);
}

void tconf_fprintf(const TCONF *tconf, FILE *stream)
{
	fprintf(stream, "EG_SIZE\tPCS_SIZE\n");
	fprintf(stream, "%4zd\t%4zd\n", tconf->eg_size, tconf->pcs_size);
	fprintf(stream, "NVEL\n");
	for(size_t g=0; g<tconf->eg_size; ++g)
		fprintf(stream, "%g\t", tconf->nvel[g]);
	fprintf(stream, "\n");
	fprintf(stream, "LAMBDAS\n");
	for(size_t p=0; p<tconf->pcs_size; ++p)
		fprintf(stream, "%g\t", tconf->lambdas[p]);
	fprintf(stream, "\n");
	fprintf(stream, "BETAS\n");
	for(size_t p=0; p<tconf->pcs_size; ++p)
		fprintf(stream, "%g\t", tconf->betas[p]);
	fprintf(stream, "\n");
	fprintf(stream, "TAU\tSTEPS\n");
	fprintf(stream, "%g\t%d\n", tconf->tau, tconf->steps);
}
