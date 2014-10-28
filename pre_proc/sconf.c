#include"sconf.h"
#include<stdlib.h>

SCONF *sconf_create(const INPUT *input)
{
	SCONF *sconf = malloc(sizeof(SCONF));
	sconf->eg_size = input->eg_size;
	size_t xm_span_size = input->xm_span_size;
	size_t ym_span_size = input->ym_span_size;
	size_t zm_span_size = input->zm_span_size;
	sconf->xm_span_size = xm_span_size;
	sconf->ym_span_size = ym_span_size;
	sconf->zm_span_size = zm_span_size;
	
	sconf->xm_mesh_size = 0;
	sconf->xspan_len = input->xspan_len;
	sconf->xspan_subdiv = input->xspan_subdiv;
	for(size_t i=0; i<xm_span_size; ++i)
		sconf->xm_mesh_size += sconf->xspan_subdiv[i];
	sconf->ym_mesh_size = 0;
	sconf->yspan_len = input->yspan_len;
	sconf->yspan_subdiv = input->yspan_subdiv;
	for(size_t j=0; j<ym_span_size; ++j)
		sconf->ym_mesh_size += sconf->yspan_subdiv[j];
	sconf->zm_mesh_size = 0;
	sconf->zspan_len = input->zspan_len;
	sconf->zspan_subdiv = input->zspan_subdiv;
	for(size_t k=0; k<zm_span_size; ++k)
		sconf->zm_mesh_size += sconf->zspan_subdiv[k];
	sconf->rt_mesh_size = 0;
	for(size_t i=0; i<xm_span_size; ++i)
		for(size_t j=0; j<ym_span_size; ++j)
			for(size_t k=0; k<zm_span_size; ++k){
				if(input->mtrl_set[i][j][k] < 0) continue;
				sconf->rt_mesh_size += input->xspan_subdiv[i]
					 	     * input->yspan_subdiv[j] 
					 	     * input->zspan_subdiv[k];
			}
	sconf->xl_bdy = input->xl_bdy;
	sconf->xr_bdy = input->xr_bdy;
	sconf->yl_bdy = input->yl_bdy;
	sconf->yr_bdy = input->yr_bdy;
	sconf->zl_bdy = input->zl_bdy;
	sconf->zr_bdy = input->zr_bdy;
	
	sconf->mtrl_set = input->mtrl_set;
	sconf->mtrllib = input->mtrllib;
	return sconf;
}

void sconf_free(SCONF *sconf)
{
	free(sconf);
}

void sconf_fprintf(const SCONF *sconf, FILE *stream)
{
	fprintf(stream, "EG_SIZE:\n");
	fprintf(stream, "%4zd\n", sconf->eg_size);
	fprintf(stream, "XM_SPAN_SIZE\tYM_SPAN_SIZE\tZM_SPAN_SIZE:\n");
	fprintf(stream, "%4zd\t%4zd\t%4zd\n", sconf->xm_span_size, sconf->ym_span_size, sconf->zm_span_size);
	fprintf(stream, "XM_MESH_SIZE\tYM_MESH_SIZE\tZM_MESH_SIZE\tRT_MESH_SIZE:\n");
	fprintf(stream, "%4zd\t%4zd\t%4zd\t%4zd\n", sconf->xm_mesh_size, sconf->ym_mesh_size, sconf->zm_mesh_size, sconf->rt_mesh_size);
	fprintf(stream, "XSPAN_LEN:\n");
	for(size_t i=0; i<sconf->xm_span_size; ++i)
		fprintf(stream, "%4g\t", sconf->xspan_len[i]);
	fprintf(stream, "\n");
	fprintf(stream, "YSPAN_LEN:\n");
	for(size_t i=0; i<sconf->ym_span_size; ++i)
		fprintf(stream, "%4g\t", sconf->yspan_len[i]);
	fprintf(stream, "\n");
	fprintf(stream, "ZSPAN_LEN:\n");
	for(size_t i=0; i<sconf->zm_span_size; ++i)
		fprintf(stream, "%4g\t", sconf->zspan_len[i]);
	fprintf(stream, "\n");
	fprintf(stream, "XSPAN_SUBDIV:\n");
	for(size_t i=0; i<sconf->xm_span_size; ++i)
		fprintf(stream, "%4zd\t", sconf->xspan_subdiv[i]);
	fprintf(stream, "\n");
	fprintf(stream, "YSPAN_SUBDIV:\n");
	for(size_t i=0; i<sconf->ym_span_size; ++i)
		fprintf(stream, "%4zd\t", sconf->yspan_subdiv[i]);
	fprintf(stream, "\n");
	fprintf(stream, "ZSPAN_SUBDIV:\n");
	for(size_t i=0; i<sconf->zm_span_size; ++i)
		fprintf(stream, "%4zd\t", sconf->zspan_subdiv[i]);
	fprintf(stream, "\n");
	fprintf(stream, "XL_BDY\tXR_BDY:\n");
	fprintf(stream, "%4d\t%4d\n", sconf->xl_bdy, sconf->xr_bdy);
	fprintf(stream, "YL_BDY\tYR_BDY:\n");
	fprintf(stream, "%4d\t%4d\n", sconf->yl_bdy, sconf->yr_bdy);
	fprintf(stream, "ZL_BDY\tZR_BDY:\n");
	fprintf(stream, "%4d\t%4d\n", sconf->zl_bdy, sconf->zr_bdy);
	fprintf(stream, "\n");
	fprintf(stream, "MATERIAL SETTING:\n");
	for(size_t k=0; k<sconf->zm_span_size; ++k){
		fprintf(stream, "z = %zd\n", k);
		for(size_t i=0; i<sconf->xm_span_size; ++i){
			for(size_t j=0; j<sconf->ym_span_size; ++j)
				fprintf(stream,"%4d\t", sconf->mtrl_set[i][j][k]);
			fprintf(stream, "\n");
		}
	}
	fprintf(stream, "\n");
	mtrllib_fprintf(sconf->mtrllib, stream);
	
}

MBLOCK sconf_get_mblock(const SCONF *sconf, size_t xspan, size_t yspan, size_t zspan)
{
	size_t xm_span_size = sconf->xm_span_size;
	size_t ym_span_size = sconf->ym_span_size;
	size_t zm_span_size = sconf->zm_span_size;
	#ifdef DEBUG
	if(xspan >= xm_span_size|| yspan >= ym_span_size || zspan >= zm_span_size){
		fprintf(stderr, "Span index out of range.\n");
		exit(-1);
	}
	#endif
	size_t *xspan_subdiv = sconf->xspan_subdiv;
	size_t *yspan_subdiv = sconf->yspan_subdiv;
	size_t *zspan_subdiv = sconf->zspan_subdiv; 
	//int ***mtrl_set = sconf->mtrl_set;
	//#ifdef DEBUG
	//if(mtrl_set[xspan][yspan][zspan] < 0){
	//	fprintf(stderr, "Specified span has no material filled in.\n");
	//	exit(-1);
	//}
	//#endif
	MBLOCK mblk;
	mblk.start_x = 0;
	for(size_t i=0; i<xspan; ++i)
		mblk.start_x += xspan_subdiv[i];
	mblk.end_x = mblk.start_x + xspan_subdiv[xspan] - 1;
	mblk.start_y = 0;
	for(size_t j=0; j<yspan; ++j)
		mblk.start_y += yspan_subdiv[j];
	mblk.end_y = mblk.start_y + yspan_subdiv[yspan] - 1;
	mblk.start_z = 0;
	for(size_t k=0; k<zspan; ++k)
		mblk.start_z += zspan_subdiv[k];
	mblk.end_z = mblk.start_z + zspan_subdiv[zspan] - 1;
	return mblk;
}
