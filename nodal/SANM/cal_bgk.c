inline int delta_func(size_t g, size_t from_g)
{
	return g == from_g? 1:0;
}


inline double cal_bgk(size_t g, size_t from_g, double *Dgk, double *srgk, double *chigk, double **ssgk, double *vsfgk, double keff)
{
	return srgk[g]*delta_func(g,from_g) - ssgk[g][from_g] - chigk[g]*vsfgk[from_g] / keff;
}
