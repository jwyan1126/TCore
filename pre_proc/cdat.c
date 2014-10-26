#include"cdat.h"
#include<stdlib.h>

CDAT3 *cdat3_create(size_t xsize, size_t ysize, size_t zsize)
{
	CDAT3 *dat = malloc(sizeof(CDAT3));
	dat->xsize = xsize;
	dat->ysize = ysize;
	dat->zsize = zsize;
	dat->data = malloc(zsize * sizeof(double **));
	for(size_t k=0; k< zsize; ++k){
		dat->data[k] = malloc(ysize * sizeof(double *));
		for(size_t j=0; j< ysize; ++j)
			dat->data[k][j] = calloc(xsize, sizeof(double));
	}
	return dat;
}

void cdat3_free(CDAT3 *dat)
{
	size_t zsize = dat->zsize;
	size_t ysize = dat->ysize;
	for(size_t k=0; k< zsize; ++k){
		for(size_t j=0; j< ysize; ++j)
			free(dat->data[k][j]);
		free(dat->data[k]);
	}
	free(dat->data);
	free(dat);
}

inline double cdat3_get_val(const CDAT3 *dat, size_t i, size_t j, size_t k)
{
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	#ifdef DEBUG
	if(i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return dat->data[k][j][i];
}

inline void cdat3_set_val(CDAT3 *dat, size_t i, size_t j, size_t k, double val)
{
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	#ifdef DEBUG
	if(i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	dat->data[k][j][i] = val;
}

void cdat3_copy(CDAT3 *tar_dat, const CDAT3 *src_dat)
{
	size_t xsize = tar_dat->xsize;
	size_t ysize = tar_dat->ysize;
	size_t zsize = tar_dat->zsize;
	#ifdef DEBUG
	if(xsize != src_dat->xsize ||
	   ysize != src_dat->ysize ||
	   zsize != src_dat->zsize){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	for(size_t k=0; k<zsize; ++k)
		for(size_t j=0; j<ysize; ++j)
			for(size_t i=0; i<xsize; ++i)
				tar_dat->data[k][j][i] = src_dat->data[k][j][i];
}

CDAT4 *cdat4_create(size_t gsize, size_t xsize, size_t ysize, size_t zsize)
{
	CDAT4 *dat = malloc(sizeof(CDAT4));
	dat->gsize = gsize;
	dat->xsize = xsize;
	dat->ysize = ysize;
	dat->zsize = zsize;
	dat->data = malloc(zsize * sizeof(double ***));
	for(size_t k=0; k< zsize; ++k){
		dat->data[k] = malloc(ysize * sizeof(double **));
		for(size_t j=0; j< ysize; ++j){
			dat->data[k][j] = malloc(xsize * sizeof(double *));
			for(size_t i=0; i< xsize; ++i)
				dat->data[k][j][i] = calloc(gsize, sizeof(double));
		}
	}
	return dat;
}

void cdat4_free(CDAT4 *dat)
{
	size_t zsize = dat->zsize;
	size_t ysize = dat->ysize;
	size_t xsize = dat->xsize;
	for(size_t k=0; k< zsize; ++k){
		for(size_t j=0; j< ysize; ++j){
			for(size_t i=0; i< xsize; ++i)
				free(dat->data[k][j][i]);
			free(dat->data[k][j]);
		}
		free(dat->data[k]);
	}
	free(dat->data);
	free(dat);
}

double cdat4_get_val(const CDAT4 *dat, size_t g, size_t i, size_t j, size_t k)
{
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	#ifdef DEBUG
	if(g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return dat->data[k][j][i][g];
}

void cdat4_set_val(CDAT4 *dat, size_t g, size_t i, size_t j, size_t k, double val)
{
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	#ifdef DEBUG
	if(g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	dat->data[k][j][i][g] = val;
}

void cdat4_copy(CDAT4 *tar_dat, const CDAT4 *src_dat)
{
	size_t gsize = tar_dat->gsize;
	size_t xsize = tar_dat->xsize;
	size_t ysize = tar_dat->ysize;
	size_t zsize = tar_dat->zsize;
	#ifdef DEBUG
	if(gsize != src_dat->gsize ||
	   xsize != src_dat->xsize ||
	   ysize != src_dat->ysize ||
	   zsize != src_dat->zsize){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	for(size_t k=0; k<zsize; ++k)
		for(size_t j=0; j<ysize; ++j)
			for(size_t i=0; i<xsize; ++i)
				for(size_t g=0; g<gsize; ++g)
					tar_dat->data[k][j][i][g] = src_dat->data[k][j][i][g];
}

CDAT5 *cdat5_create(size_t gsize, size_t xsize, size_t ysize, size_t zsize)
{
	CDAT5 *dat = malloc(sizeof(CDAT5));
	dat->gsize = gsize;
	dat->xsize = xsize;
	dat->ysize = ysize;
	dat->zsize = zsize;
	dat->data = malloc(zsize * sizeof(double ****));
	for(size_t k=0; k< zsize; ++k){
		dat->data[k] = malloc(ysize * sizeof(double ***));
		for(size_t j=0; j< ysize; ++j){
			dat->data[k][j] = malloc(xsize * sizeof(double **));
			for(size_t i=0; i< xsize; ++i){
				dat->data[k][j][i] = malloc(gsize * sizeof(double *));
				for(size_t g=0; g< gsize; ++g)
					dat->data[k][j][i][g] = calloc(gsize, sizeof(double));
			}
		}
	}
	return dat;
}

void cdat5_free(CDAT5 *dat)
{
	size_t gsize = dat->gsize;
	size_t zsize = dat->zsize;
	size_t ysize = dat->ysize;
	size_t xsize = dat->xsize;
	for(size_t k=0; k< zsize; ++k){
		for(size_t j=0; j< ysize; ++j){
			for(size_t i=0; i< xsize; ++i){
				for(size_t g=0; g< gsize; ++g)
					free(dat->data[k][j][i][g]);
				free(dat->data[k][j][i]);
			}
			free(dat->data[k][j]);
		}
		free(dat->data[k]);
	}
	free(dat->data);
	free(dat);
}

double cdat5_get_val(const CDAT5 *dat, size_t g, size_t from_g, size_t i, size_t j, size_t k)
{
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	#ifdef DEBUG
	if(g >= gsize || from_g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return dat->data[k][j][i][g][from_g];
}

void cdat5_set_val(CDAT5 *dat, size_t g, size_t from_g, size_t i, size_t j, size_t k, double val)
{
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	#ifdef DEBUG
	if(g >= gsize || from_g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	dat->data[k][j][i][g][from_g] = val;
}

void cdat5_copy(CDAT5 *tar_dat, const CDAT5 *src_dat)
{
	size_t gsize = tar_dat->gsize;
	size_t xsize = tar_dat->xsize;
	size_t ysize = tar_dat->ysize;
	size_t zsize = tar_dat->zsize;
	#ifdef DEBUG
	if(gsize != src_dat->gsize ||
	   xsize != src_dat->xsize ||
	   ysize != src_dat->ysize ||
	   zsize != src_dat->zsize){
		fprintf(stderr, "Size incompatible.\n");
		exit(-1);
	}
	#endif
	for(size_t k=0; k<zsize; ++k)
		for(size_t j=0; j<ysize; ++j)
			for(size_t i=0; i<xsize; ++i)
				for(size_t g=0; g<gsize; ++g)
					for(size_t from_g=0; from_g<gsize; ++from_g)
						tar_dat->data[k][j][i][g][from_g] = src_dat->data[k][j][i][g][from_g];
	
}
