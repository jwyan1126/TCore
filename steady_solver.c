#include"steady_solver.h"
#include"ssol.h"
#include"algebra/mat.h"
#include"algebra/vec.h"
#include"algebra/eigen.h"
#include<stdlib.h>
#include"nodal/leak.h"
#include"nodal/jcur.h"
#include"nodal/tnsol.h"
#include"nodal/SANM/sanm.h"

void steady_state_cal(SSOL *ssol, SCONF *sconf, MAPPER *mapper, const MESH *mesh)
{
	size_t eg_size = sconf->eg_size;
	size_t rt_size = sconf->rt_mesh_size;
	MAT *M = mat_create(eg_size * rt_size);
	MAT *S = mat_create(eg_size * rt_size);
	MAT *F = mat_create(eg_size * rt_size);
	EDAT4 *DFDM = edat4_create(mapper);
	EDAT4 *DNOD = edat4_create(mapper);
	EDAT4 *Jn = edat4_create(mapper);
	LEAK *leak = leak_create(mapper);
	VEC *phi = vec_ref_create(eg_size * rt_size, ssol->flux->data);
	vec_normalize(phi);
	VEC *phi_last = vec_create(eg_size * rt_size);
	
	cal_DFDM(DFDM, mesh);
	cal_s(S, sconf,mapper, mesh);
	cal_f(F, sconf, mapper, mesh);
	cal_m(M, DFDM, DNOD, mapper, mesh, NULL);
	mat_adds(M, S, -1.0);
	double k;
	gspow_iter(&k, phi, M, F, 0.0, 2048);
	ssol->keff = 1.0 / k;
	vec_free(phi_last);
	leak_free(leak);
	edat4_free(Jn);
	edat4_free(DNOD);
	edat4_free(DFDM);
	vec_ref_free(phi);
	mat_free(F);
	mat_free(S);
	mat_free(M);
}
