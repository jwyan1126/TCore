SCONF *sconf_create(INPUT *input)
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
		sconf->xm_mesh_size += xspan_subdiv[i];
	sconf->ym_mesh_size = 0;
	sconf->yspan_len = input->yspan_len;
	sconf->yspan_subdiv = input->yspan_subdiv;
	for(size_t j=0; j<ym_span_size; ++j)
		sconf->ym_mesh_size += yspan_subdiv[j];
	sconf->zm_mesh_size = 0;
	sconf->zspan_len = input->zspan_len;
	sconf->zspan_subdiv = input->zspan_subdiv;
	for(size_t k=0; k<zm_span_size; ++k)
		sconf->zm_mesh_size += zspan_subdiv[k];
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
