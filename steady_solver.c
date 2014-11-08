#include"steady_solver.h"
#include"algebra/mat.h"
#include"algebra/vec.h"
#include"algebra/eigen.h"
#include<stdlib.h>
#include"nodal/leak.h"
#include"nodal/jcur.h"
#include"nodal/tnsol.h"
#include"nodal/SANM/sanm.h"

double steady_state_cal(FLUX *flux, SCONF *sconf, MAPPER *mapper, const MESH *mesh)
{
	size_t eg_size = sconf->eg_size;
	size_t rt_size = sconf->rt_mesh_size;
	MAT *M = mat_create(eg_size * rt_size);
	MAT *S = mat_create(eg_size * rt_size);
	MAT *F = mat_create(eg_size * rt_size);
	EDAT4 *DFDM = edat4_create(mapper);
	EDAT4 *DNOD = edat4_create(mapper);
	VEC *phi = vec_ref_create(eg_size * rt_size, flux->data);
	cal_DFDM(DFDM, mesh);
	cal_s(S, mapper, mesh);
	cal_f(F, mapper, mesh, NULL);
	cal_m(M, DFDM, DNOD, mapper, mesh, NULL);
	mat_adds(M, S, -1.0);
	double k;
	gspow_iter(&k, phi, M, F, 0.0, 1024);
	edat4_free(DNOD);
	edat4_free(DFDM);
	vec_ref_free(phi);
	mat_free(F);
	mat_free(S);
	mat_free(M);
	return 1.0/k;
}
