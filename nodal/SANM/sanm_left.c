#include"sanm.h"

void sanm_left(TNSOL *tn)
{
	size_t eg_size = tn->eg_size;
	double duj = tn->duj;
	double *Dgj = tn->Dgj;
	double *vsfgj = tn->vsfgj;
	double *phigj = tn->phigj;
	double *srgj = tn->srgj;
	double **ssgj = tn->ssgj;
	double *chigj = tn->chigj;
	double *lgj0 = tn->lgj0;
	double *lgj1 = tn->lgj1;
	double *lgj2 = tn->lgj2;
	double keff = tn->keff;
	double *J = tn->J;
	double *dgxj = malloc(eg_size * sizeof(double));
	double *agj1 = malloc(eg_size * sizeof(double));
	double *agj2 = malloc(eg_size * sizeof(double));
	double *agj3 = malloc(eg_size * sizeof(double));
	double *agj5 = malloc(eg_size * sizeof(double));
	double *alpha_gxj = malloc(eg_size * sizeof(double));
	double *mgxj0 = malloc(eg_size * sizeof(double));
	double *mgxj1 = malloc(eg_size * sizeof(double));
	double *mgxj2 = malloc(eg_size * sizeof(double));
	double *Agxj = malloc(eg_size * sizeof(double));
	double *Cgxj = malloc(eg_size * sizeof(double));
	double *Egxj = malloc(eg_size * sizeof(double));
	double *Fgxj = malloc(eg_size * sizeof(double));
	double *Ggxj = malloc(eg_size * sizeof(double));
	double *Hgxj = malloc(eg_size * sizeof(double));
	
	for(size_t g=0; g<eg_size; ++g){
		dgxj[g] = 2.0*Dgj[g] / duj;
		alpha_gxj[g] = sqrt(srgj[g]/Dgj[g])*duj/2.0;
    		mgxj0[g] = sinh(alpha_gxj[g])/alpha_gxj[g];
    		mgxj1[g] = 3.0*(alpha_gxj[g]*cosh(alpha_gxj[g])-sinh(alpha_gxj[g]))/(alpha_gxj[g]*alpha_gxj[g]);
    		mgxj2[g] = 5.0*(alpha_gxj[g]*alpha_gxj[g]*sinh(alpha_gxj[g])-3.0*alpha_gxj[g]*cosh(alpha_gxj[g])+3.0*sinh(alpha_gxj[g]))/(alpha_gxj[g]*alpha_gxj[g]*alpha_gxj[g]);
   		Agxj[g] = (sinh(alpha_gxj[g])-mgxj1[g])/(srgj[g]*mgxj1[g]);
    		Cgxj[g] = (cosh(alpha_gxj[g])-mgxj0[g]-mgxj2[g])/(srgj[g]*mgxj2[g]);
    		Egxj[g] = (mgxj0[g]/mgxj2[g])-3.0/(alpha_gxj[g]*alpha_gxj[g]);
    		Fgxj[g] = (alpha_gxj[g]*cosh(alpha_gxj[g])-mgxj1[g])/(srgj[g]*mgxj1[g]);
    		Ggxj[g] = (alpha_gxj[g]*sinh(alpha_gxj[g])-3.0*mgxj2[g])/(cosh(alpha_gxj[g])-mgxj0[g]-mgxj2[g]);
    		Hgxj[g] = (alpha_gxj[g]*cosh(alpha_gxj[g])-mgxj1[g])/(sinh(alpha_gxj[g])-mgxj1[g]);
	}

	// agi2
	MAT *A = mat_create(eg_size);
	for(size_t g=0; g<eg_size; ++g)
		for(size_t from_g=0; from_g<eg_size; ++from_g){
			double val = 12.0*Dgi[g]*delta_func(g,from_g)/(dui*dui) + Egxi[g]*cal_bgk(g,from_g,Dgi,srgi,chigi,ssgi,vsfgi,keff);
			mat_set(A, g, from_g, val);
		}
	VEC *b = vec_create(eg_size);
	for(size_t g=0; g<eg_size; ++g){
		double val = 0.0;
		for(size_t from_g=0; from_g<eg_size; ++from_g)
			val += cal_bgk(g,from_g,Dgi,srgi,chigi,ssgi,vsfgi,keff)*phigi[from_g];
		val += lgi0[g] - Egxi[g]*lgi2[g];
		vec_set(b, g, val);
	}
	VEC *x = vec_create(eg_size);
	gauss_seidel(x, A, b, INT_MAX);
	for(size_t g=0; g<eg_size; ++g)
		agi2[g] = vec_get(x, g);
	vec_free(x);
	vec_free(b);
	mat_free(A);
	// agi4
	for(size_t g=0; g<eg_size; ++g){
		double agi4[g] = 0.0;
		for(size_t from_g=0; from_g<eg_size; ++from_g)
			agi4[g] = cal_bgk(g,from_g,Dgi,srgi,chigi,ssgi,vsfgi,keff)*agi2[from_g];
		agi4[g] += lgi2[g];
		agi4[g] *= Cgxi[g];
	}
	// agi1
	
}
