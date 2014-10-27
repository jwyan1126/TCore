#include"edat.h"
#include<stdlib.h>

EDAT4 *edat4_create(size_t gsize, size_t xsize, size_t ysize, size_t zsize)
{
	EDAT4 *dat = malloc(sizeof(EDAT4));
	dat->gsize = gsize;
	dat->xsize = xsize;
	dat->ysize = ysize;
	dat->zsize = zsize;
	dat->xdata = malloc(zsize * sizeof(double ***));
	for(size_t k=0; k< zsize; ++k){
		dat->xdata[k] = malloc(ysize * sizeof(double **));
		for(size_t j=0; j< ysize; ++j){
			dat->xdata[k][j] = malloc((xsize+1) * sizeof(double *));
			for(size_t i=0; i< xsize+1; ++i)
				dat->xdata[k][j][i] = calloc(gsize, sizeof(double));
		}
	}
	dat->ydata = malloc(xsize * sizeof(double ***));
	for(size_t i=0; i< xsize; ++i){
		dat->ydata[i] = malloc(zsize * sizeof(double **));
		for(size_t k=0; k< zsize; ++k){
			dat->ydata[i][k] = malloc((ysize+1) * sizeof(double *));
			for(size_t j=0; j< ysize+1; ++j)
				dat->ydata[i][k][j] = calloc(gsize, sizeof(double));
		}
	}
	dat->zdata = malloc(ysize * sizeof(double ***));
	for(size_t j=0; j< ysize; ++j){
		dat->zdata[j] = malloc(xsize * sizeof(double **));
		for(size_t i=0; i< xsize; ++i){
			dat->zdata[j][i] = malloc((zsize+1) * sizeof(double *));
			for(size_t k=0; k< zsize+1; ++k)
				dat->zdata[j][i][k] = calloc(gsize, sizeof(double));
		}
	}
	return dat;
}

void edat4_free(EDAT4 *dat)
{
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	for(size_t k=0; k< zsize; ++k){
		for(size_t j=0; j< ysize; ++j){
			for(size_t i=0; i< xsize+1; ++i)
				free(dat->xdata[k][j][i]);
			free(dat->xdata[k][j]);
		}
		free(dat->xdata[k]);
	}
	free(dat->xdata);
	for(size_t i=0; i< xsize; ++i){
		for(size_t k=0; k< zsize; ++k){
			for(size_t j=0; j< ysize+1; ++j)
				free(dat->ydata[i][k][j]);
			free(dat->ydata[i][k]);
		}
		free(dat->ydata[i]);
	}
	free(dat->ydata);
	for(size_t j=0; j< ysize; ++j){
		for(size_t i=0; i< xsize; ++i){
			for(size_t k=0; k< zsize+1; ++k)
				free(dat->zdata[j][i][k]);
			free(dat->zdata[j][i]);
		}
		free(dat->zdata[j]);
	}
	free(dat->zdata);
	free(dat);
}

inline double edat4_get_xlval(const EDAT4 *dat, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	if(g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return dat->xdata[k][j][i][g];
}
inline double edat4_get_xrval(const EDAT4 *dat, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	if(g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return dat->xdata[k][j][i+1][g];
}

inline double edat4_get_ylval(const EDAT4 *dat, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	if(g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return dat->ydata[i][k][j][g];
}
inline double edat4_get_yrval(const EDAT4 *dat, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	if(g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return dat->ydata[i][k][j+1][g];
}
inline double edat4_get_zlval(const EDAT4 *dat, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	if(g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return dat->zdata[j][i][k][g];
}
inline double edat4_get_zrval(const EDAT4 *dat, size_t g, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	if(g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	return dat->zdata[j][i][k+1][g];
}

inline void edat4_set_xlval(EDAT4 *dat, size_t g, size_t i, size_t j, size_t k, double val)
{
	#ifdef DEBUG
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	if(g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	dat->xdata[k][j][i][g] = val;
}
inline void edat4_set_xrval(EDAT4 *dat, size_t g, size_t i, size_t j, size_t k, double val)
{
	#ifdef DEBUG
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	if(g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	dat->xdata[k][j][i+1][g] = val;
}
inline void edat4_set_ylval(EDAT4 *dat, size_t g, size_t i, size_t j, size_t k, double val)
{
	#ifdef DEBUG
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	if(g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	dat->ydata[i][k][j][g] = val;
}
inline void edat4_set_yrval(EDAT4 *dat, size_t g, size_t i, size_t j, size_t k, double val)
{
	#ifdef DEBUG
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	if(g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	dat->ydata[i][k][j+1][g] = val;
}
inline void edat4_set_zlval(EDAT4 *dat, size_t g, size_t i, size_t j, size_t k, double val)
{
	#ifdef DEBUG
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	if(g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	dat->zdata[j][i][k][g] = val;
}
inline void edat4_set_zrval(EDAT4 *dat, size_t g, size_t i, size_t j, size_t k, double val)
{
	#ifdef DEBUG
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	if(g >= gsize || i >= xsize || j >= ysize || k >= zsize){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	dat->zdata[j][i][k+1][g] = val;
}

void edat4_copy(EDAT4 *tar_dat, const EDAT4 *src_dat)
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
			for(size_t i=0; i<xsize+1; ++i)
				for(size_t g=0; g<gsize; ++g)
					tar_dat->xdata[k][j][i][g] = src_dat->xdata[k][j][i][g];
	for(size_t i=0; i<xsize; ++i)
		for(size_t k=0; k<zsize; ++k)
			for(size_t j=0; j<ysize+1; ++j)
				for(size_t g=0; g<gsize; ++g)
					tar_dat->ydata[i][k][j][g] = src_dat->ydata[i][k][j][g];
	for(size_t j=0; j<ysize; ++j)
		for(size_t i=0; i<xsize; ++i)
			for(size_t k=0; k<zsize+1; ++k)
				for(size_t g=0; g<gsize; ++g)
					tar_dat->zdata[j][i][k][g] = src_dat->zdata[j][i][k][g];
}

void edat4_fprintf(EDAT4 *dat, size_t g, size_t i, size_t j, size_t k, FILE *stream)
{
	fprintf(stream, "EG=%4zd\tX=%4zd\tY=%4zd\tZ=%4zd\n",g,i,j,k);
	fprintf(stream, "XL\t\tXR\t\tYL\t\tYR\t\tZL\t\tZR\n");
	fprintf(stream, "%4g\t%4g\t%4g\t%4g\t%4g\t%4g\n", 
			edat4_get_xlval(dat,g,i,j,k),
			edat4_get_xrval(dat,g,i,j,k),
			edat4_get_ylval(dat,g,i,j,k),
			edat4_get_yrval(dat,g,i,j,k),
			edat4_get_zlval(dat,g,i,j,k),
			edat4_get_zrval(dat,g,i,j,k));
}

void edat4_set_rand(EDAT4 *dat)
{
	size_t gsize = dat->gsize;
	size_t xsize = dat->xsize;
	size_t ysize = dat->ysize;
	size_t zsize = dat->zsize;
	for(size_t k=0; k<zsize; ++k)
		for(size_t j=0; j<ysize; ++j)
			for(size_t i=0; i<xsize+1; ++i)
				for(size_t g=0; g<gsize; ++g)
					dat->xdata[k][j][i][g] = rand();
	for(size_t i=0; i<xsize; ++i)
		for(size_t k=0; k<zsize; ++k)
			for(size_t j=0; j<ysize+1; ++j)
				for(size_t g=0; g<gsize; ++g)
					dat->ydata[i][k][j][g] = rand();
	for(size_t j=0; j<ysize; ++j)
		for(size_t i=0; i<xsize; ++i)
			for(size_t k=0; k<zsize+1; ++k)
				for(size_t g=0; g<gsize; ++g)
					dat->zdata[j][i][k][g] = rand();
}
