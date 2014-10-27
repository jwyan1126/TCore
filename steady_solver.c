#include"solver.h"
#include"solution.h"
#include"algebra/mat.h"
#include"algebra/vec.h"

void steady_solver(SSOL *ssol, const SCONF *sconf, const MESH *mesh)
{
	size_t eg_size = sconf->eg_size;
	size_t xm_size = sconf->xm_mesh_size;
	size_t ym_size = sconf->ym_mesh_size;
	size_t zm_size = sconf->zm_mesh_size;
	size_t rt_size = sconf->rt_mesh_size;
	MAT *M = mat_create(eg_size * rt_size);
	MAT *S = mat_create(eg_size * rt_size);
	MAT *F = mat_create(eg_size * rt_size);
	EDAT4 *DFDM = edat4_create(eg_size, xm_size, ym_size, zm_size);
	EDAT4 *DNOD = edat4_create(eg_size, xm_size, ym_size, zm_size);
	cal_DFDM(DFDM, sconf, mesh);
	cal_M(M, DFDM, DNOD);
	cal_S(S, sconf, mesh);
	mat_adds(M, S, -1.0);
	cal_F(F, sconf, mesh);
	VEC *phi = vec_ref_create(eg_size * rt_size, ssol->flux);
	if(!gspow_iter(&ssol->keff, phi, M, F, 0.0, 128))
		printf("Steady cal. converged.\n");
	else
		printf("Steady cal. NOT converged.\n");
	edat4_free(DNOD);
	edat4_free(DFDM);
	vec_free(phi);
	mat_free(F);
	mat_free(S);
	mat_free(M);
}

void cal_DFDM(EDAT *DFDM, const *sconf, const *mesh)
{
	size_t eg_size = sconf->eg_size;
	size_t xm_size = sconf->xm_mesh_size;
	size_t ym_size = sconf->ym_mesh_size;
	size_t zm_size = sconf->zm_mesh_size;
	int ***bdy_checker = mesh->bdy_checker;
	double xl_bdy = sconf->xl_bdy;
	double xr_bdy = sconf->xr_bdy;
	double yl_bdy = sconf->yl_bdy;
	double yr_bdy = sconf->yr_bdy;
	double zl_bdy = sconf->zl_bdy;
	double zr_bdy = sconf->zr_bdy;
	for(size_t k=0; k<zm_size; ++k)
		for(size_t j=0; j<ym_size; ++j)
			for(size_t i=0; i<xm_size; ++i){
				//xdata[k][j][i]
				if(bdy_checker[k][j][i] & 0b00000001) continue;
				if()//.................
			}
	for(size_t i=0; i<xm_size; ++i)
		for(size_t k=0; k<zm_size; ++k)
			for(size_t j=0; j<ym_size; ++j){
				//ydata[i][k][j]
			}
	for(size_t j=0; j<ym_size; ++j)
		for(size_t i=0; i<xm_size; ++i)
			for(size_t k=0; k<zm_size; ++k){
				//zdata[j][i][k]
			}
}

void cal_M(MAT *M, )
{
	
}

void cal_S(MAT *S, const *sconf, const *mesh)
{
}

void cal_F(MAT *F, const *sconf, const *mesh)
{

}
