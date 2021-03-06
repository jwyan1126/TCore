#include"mapper.h"
#include<stdlib.h>
#include<stdio.h>
#include"checker.h"

size_t inwhichspan(size_t arr[], size_t span_len, int i);

MAPPER *mapper_create(SCONF *sconf)
{
	size_t eg_size = sconf->eg_size;
	size_t xm_size = sconf->xm_mesh_size;
	size_t ym_size = sconf->ym_mesh_size;
	size_t zm_size = sconf->zm_mesh_size;
	size_t rt_size = sconf->rt_mesh_size;
	MAPPER *mapper = malloc(sizeof(MAPPER));
	mapper->eg_size = eg_size;
	mapper->xm_size = xm_size;
	mapper->ym_size = ym_size;
	mapper->zm_size = zm_size;
	mapper->rt_size = rt_size;
	mapper->one2three = malloc(rt_size * sizeof(XYZ_IDX));
	mapper->three2one = malloc(zm_size * sizeof(size_t **));
	mapper->cchecker = malloc(zm_size * sizeof(int **));
	for(size_t k=0; k<zm_size; ++k){
		mapper->cchecker[k] = malloc(ym_size * sizeof(int *));
		mapper->three2one[k] = malloc(ym_size * sizeof(size_t *));
		for(size_t j=0; j<ym_size; ++j){
			mapper->cchecker[k][j] = calloc(xm_size, sizeof(int));
			mapper->three2one[k][j] = calloc(xm_size, sizeof(size_t));
		}
	}
	mapper->xchecker = malloc(zm_size * sizeof(int **));
	for(size_t k=0; k< zm_size; ++k){
		mapper->xchecker[k] = malloc(ym_size * sizeof(int *));
		for(size_t j=0; j<ym_size; ++j){
			mapper->xchecker[k][j] = calloc(xm_size+1, sizeof(int));
		}
	}
	mapper->ychecker = malloc(xm_size * sizeof(int **));
	for(size_t i=0; i<xm_size; ++i){
		mapper->ychecker[i] = malloc(zm_size * sizeof(int *));
		for(size_t k=0; k<zm_size; ++k){
			mapper->ychecker[i][k] = calloc(ym_size+1, sizeof(int));
		}
	}
	mapper->zchecker = malloc(ym_size * sizeof(int **));
	for(size_t j=0; j< ym_size; ++j){
		mapper->zchecker[j] = malloc(xm_size * sizeof(int *));
		for(size_t i=0; i< xm_size; ++i){
			mapper->zchecker[j][i] = calloc(zm_size+1, sizeof(int));
		}
	}
	// traversal
	size_t ac = 0;
	for(size_t k=0; k< zm_size; ++k)
		for(size_t j=0; j< ym_size; ++j)
			for(size_t i=0; i< xm_size; ++i){
				size_t xspan = inwhichspan(sconf->xspan_subdiv, sconf->xm_mesh_size, i);
				size_t yspan = inwhichspan(sconf->yspan_subdiv, sconf->ym_mesh_size, j);
				size_t zspan = inwhichspan(sconf->zspan_subdiv, sconf->zm_mesh_size, k);
				int mtrl_id = sconf->mtrl_set[xspan][yspan][zspan];
				if(mtrl_id < 0)
					continue;
				mapper->three2one[k][j][i] = ac;
				XYZ_IDX xyz;
				xyz.xi = i; xyz.yi = j; xyz.zi = k;
				mapper->one2three[ac] = xyz;
				ac++;
			}
	bdy_check(mapper->cchecker, sconf);
	bdy_check_edge(mapper->xchecker, mapper->ychecker, mapper->zchecker, sconf);
	return mapper;
}

void mapper_free(MAPPER *mapper)
{
	size_t xm_size = mapper->xm_size;
	size_t zm_size = mapper->zm_size;
	size_t ym_size = mapper->ym_size;
	for(size_t k=0; k<mapper->zm_size; ++k){
		for(size_t j=0; j<mapper->ym_size; ++j){
			free(mapper->three2one[k][j]);
			free(mapper->cchecker[k][j]);
		}
		free(mapper->three2one[k]);
		free(mapper->cchecker[k]);
	}
	for(size_t k=0; k< zm_size; ++k){
		for(size_t j=0; j< ym_size; ++j){
			free(mapper->xchecker[k][j]);
		}
		free(mapper->xchecker[k]);
	}
	free(mapper->xchecker);
	for(size_t i=0; i< xm_size; ++i){
		for(size_t k=0; k< zm_size; ++k){
			free(mapper->ychecker[i][k]);
		}
		free(mapper->ychecker[i]);
	}
	free(mapper->ychecker);
	for(size_t j=0; j< ym_size; ++j){
		for(size_t i=0; i< xm_size; ++i){
			free(mapper->zchecker[j][i]);
		}
		free(mapper->zchecker[j]);
	}
	free(mapper->zchecker);
	free(mapper->three2one);
	free(mapper->one2three);
	free(mapper->cchecker);
	free(mapper);
}

inline XYZ_IDX mapper_get3Didx(const MAPPER *mapper, size_t idx1D)
{
	#ifdef DEBUG
	size_t rt_size = mapper->rt_size;
	if(idx1D >= rt_size){
		fprintf(stderr, "'idx1D' out of range.\n");
		exit(-1);
	}
	#endif
	return mapper->one2three[idx1D];
}

inline size_t mapper_get1Didx(const MAPPER *mapper, size_t i, size_t j, size_t k)
{
	#ifdef DEBUG
	size_t xm_size = mapper->xm_size;
	size_t ym_size = mapper->ym_size;
	size_t zm_size = mapper->zm_size;
	if(i >= xm_size || j >= ym_size || k >= zm_size){
		fprintf(stderr, "'idx3D' out of range.\n");
		exit(-1);
	}
	#endif
	return mapper->three2one[k][j][i];
}

void mapper_fprintf(const MAPPER *mapper, FILE *stream)
{
	fprintf(stream, "1DPOS\t3DXPOS\t3DYPOS\t3DZPOS\n");
	for(size_t i=0; i < mapper->rt_size; ++i){
		size_t x = mapper->one2three[i].xi;
		size_t y = mapper->one2three[i].yi;
		size_t z = mapper->one2three[i].zi;
		size_t r = mapper->three2one[z][y][x];
		fprintf(stream, "%zd\t%zd\t%zd\t%zd\n",r,x,y,z);
	}
}

// private func.
size_t inwhichspan(size_t arr[], size_t span_len, int i)
{
	for(size_t k=0; k<span_len; ++k){
		i -= arr[k];
		if(i < 0)
			return k;
	}
	fprintf(stderr, "Error occurs. when invoke 'inwhichspan'.\n");
	exit(-1);
}
