#ifndef INPUT_H
#define INPUT_H

typedef struct
{
	size_t eg_size;
	size_t xm_span_size;
	size_t ym_span_size;
	size_t zm_span_size;
	double *xspan_len;
	double *yspan_len;
	double *zspan_len;
	size_t *xspan_subdiv;
	size_t *yspan_subdiv;
	size_t *zspan_subdiv;

	double xl_bdy;
	double xr_bdy;
	double yl_bdy;
	double yr_bdy;
	double zl_bdy;
	double zr_bdy;

	int ***mtrl_set;
	MTRLLIB *mtrllib;
} INPUT;

INPUT *input_create(const char *path);

input_free(INPUT *input);

#endif
