#include"sanm.h"
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"../../algebra/mat.h"
#include"../../algebra/vec.h"
#include"../../algebra/ksp.h"
#include"cal_bgk.h"
#include<limits.h>

void sanm_inner(TNSOL *tn)
{
	size_t eg_size = tn->eg_size;
	double dui = tn->dui;  		double duj = tn->duj;
	double *Dgi = tn->Dgi; 		double *Dgj = tn->Dgj;
	double *vsfgi = tn->vsfgi;	double *vsfgj = tn->vsfgj;
	double *phigi = tn->phigi;	double *phigj = tn->phigj;
	double *srgi = tn->srgi;	double *srgj = tn->srgj;
	double **ssgi = tn->ssgi;	double **ssgj = tn->ssgj;
	double *chigi = tn->chigi;	double *chigj = tn->chigj;
	double *lgi0 = tn->lgi0;	double *lgj0 = tn->lgj0;
	double *lgi1 = tn->lgi1;	double *lgj1 = tn->lgj1;
	double *lgi2 = tn->lgi2;	double *lgj2 = tn->lgj2;
	double *adfgi = tn->adfgi;	double *adfgj = tn->adfgj;
	double keff = tn->keff;		double *J = tn->J;
	double *dgxi = malloc(eg_size * sizeof(double));
	double *dgxj = malloc(eg_size * sizeof(double));
	double *agi1 = malloc(eg_size * sizeof(double));
	double *agj1 = malloc(eg_size * sizeof(double));
	double *agi2 = malloc(eg_size * sizeof(double));
	double *agj2 = malloc(eg_size * sizeof(double));
	double *agi3 = malloc(eg_size * sizeof(double));
	double *agj3 = malloc(eg_size * sizeof(double));
	double *agi4 = malloc(eg_size * sizeof(double));
	double *agj4 = malloc(eg_size * sizeof(double));
	double *alpha_gxi = malloc(eg_size * sizeof(double));
	double *alpha_gxj = malloc(eg_size * sizeof(double));
	double *mgxi0 = malloc(eg_size * sizeof(double));
	double *mgxj0 = malloc(eg_size * sizeof(double));
	double *mgxi1 = malloc(eg_size * sizeof(double));
	double *mgxj1 = malloc(eg_size * sizeof(double));
	double *mgxi2 = malloc(eg_size * sizeof(double));
	double *mgxj2 = malloc(eg_size * sizeof(double));
	double *Agxi = malloc(eg_size * sizeof(double));
	double *Agxj = malloc(eg_size * sizeof(double));
	double *Cgxi = malloc(eg_size * sizeof(double));
	double *Cgxj = malloc(eg_size * sizeof(double));
	double *Egxi = malloc(eg_size * sizeof(double));
	double *Egxj = malloc(eg_size * sizeof(double));
	double *Fgxi = malloc(eg_size * sizeof(double));
	double *Fgxj = malloc(eg_size * sizeof(double));
	double *Ggxi = malloc(eg_size * sizeof(double));
	double *Ggxj = malloc(eg_size * sizeof(double));
	double *Hgxi = malloc(eg_size * sizeof(double));
	double *Hgxj = malloc(eg_size * sizeof(double));
	for(size_t g=0; g<eg_size; ++g){
		dgxi[g] = 2.0*Dgi[g] / dui;
		dgxj[g] = 2.0*Dgj[g] / duj;
		alpha_gxi[g] = sqrt(srgi[g]/Dgi[g])*dui/2.0;
		alpha_gxj[g] = sqrt(srgj[g]/Dgj[g])*duj/2.0;
    		mgxi0[g] = sinh(alpha_gxi[g])/alpha_gxi[g];
    		mgxj0[g] = sinh(alpha_gxj[g])/alpha_gxj[g];
    		mgxi1[g] = 3.0*(alpha_gxi[g]*cosh(alpha_gxi[g])-sinh(alpha_gxi[g]))/(alpha_gxi[g]*alpha_gxi[g]); 
    		mgxj1[g] = 3.0*(alpha_gxj[g]*cosh(alpha_gxj[g])-sinh(alpha_gxj[g]))/(alpha_gxj[g]*alpha_gxj[g]);
    		mgxi2[g] = 5.0*(alpha_gxi[g]*alpha_gxi[g]*sinh(alpha_gxi[g])-3.0*alpha_gxi[g]*cosh(alpha_gxi[g])+3.0*sinh(alpha_gxi[g]))/(alpha_gxi[g]*alpha_gxi[g]*alpha_gxi[g]);
    		mgxj2[g] = 5.0*(alpha_gxj[g]*alpha_gxj[g]*sinh(alpha_gxj[g])-3.0*alpha_gxj[g]*cosh(alpha_gxj[g])+3.0*sinh(alpha_gxj[g]))/(alpha_gxj[g]*alpha_gxj[g]*alpha_gxj[g]);
    		Agxi[g] = (sinh(alpha_gxi[g])-mgxi1[g])/(srgi[g]*mgxi1[g]);
   		Agxj[g] = (sinh(alpha_gxj[g])-mgxj1[g])/(srgj[g]*mgxj1[g]);
    		Cgxi[g] = (cosh(alpha_gxi[g])-mgxi0[g]-mgxi2[g])/(srgi[g]*mgxi2[g]);
    		Cgxj[g] = (cosh(alpha_gxj[g])-mgxj0[g]-mgxj2[g])/(srgj[g]*mgxj2[g]);
    		Egxi[g] = (mgxi0[g]/mgxi2[g])-3.0/(alpha_gxi[g]*alpha_gxi[g]);
    		Egxj[g] = (mgxj0[g]/mgxj2[g])-3.0/(alpha_gxj[g]*alpha_gxj[g]);
    		Fgxi[g] = (alpha_gxi[g]*cosh(alpha_gxi[g])-mgxi1[g])/(srgi[g]*mgxi1[g]);
    		Fgxj[g] = (alpha_gxj[g]*cosh(alpha_gxj[g])-mgxj1[g])/(srgj[g]*mgxj1[g]);
    		Ggxi[g] = (alpha_gxi[g]*sinh(alpha_gxi[g])-3.0*mgxi2[g])/(cosh(alpha_gxi[g])-mgxi0[g]-mgxi2[g]);
    		Ggxj[g] = (alpha_gxj[g]*sinh(alpha_gxj[g])-3.0*mgxj2[g])/(cosh(alpha_gxj[g])-mgxj0[g]-mgxj2[g]);
    		Hgxi[g] = (alpha_gxi[g]*cosh(alpha_gxi[g])-mgxi1[g])/(sinh(alpha_gxi[g])-mgxi1[g]);
    		Hgxj[g] = (alpha_gxj[g]*cosh(alpha_gxj[g])-mgxj1[g])/(sinh(alpha_gxj[g])-mgxj1[g]);
	}
	// agi2
	{
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
	LU_solve(x, A, b);
	for(size_t g=0; g<eg_size; ++g)
		agi2[g] = vec_get(x, g);
	vec_free(x);
	vec_free(b);
	mat_free(A);
	}
	
	// agj2
	{
	MAT *A = mat_create(eg_size);
	for(size_t g=0; g<eg_size; ++g)
		for(size_t from_g=0; from_g<eg_size; ++from_g){
			double val = 12.0*Dgj[g]*delta_func(g,from_g)/(duj*duj) + Egxj[g]*cal_bgk(g,from_g,Dgj,srgj,chigj,ssgj,vsfgj,keff);
			mat_set(A, g, from_g, val);
		}
	VEC *b = vec_create(eg_size);
	for(size_t g=0; g<eg_size; ++g){
		double val = 0.0;
		for(size_t from_g=0; from_g<eg_size; ++from_g)
			val += cal_bgk(g,from_g,Dgj,srgj,chigj,ssgj,vsfgj,keff)*phigj[from_g];
		val += lgj0[g] - Egxj[g]*lgj2[g];
		vec_set(b, g, val);
	}
	VEC *x = vec_create(eg_size);
	LU_solve(x, A, b);
	for(size_t g=0; g<eg_size; ++g)
		agj2[g] = vec_get(x, g);
	vec_free(x);
	vec_free(b);
	mat_free(A);
	}
	// agi4
	for(size_t g=0; g<eg_size; ++g){
		agi4[g] = 0.0;
		for(size_t from_g=0; from_g<eg_size; ++from_g)
			agi4[g] += cal_bgk(g,from_g,Dgi,srgi,chigi,ssgi,vsfgi,keff)*agi2[from_g];
		agi4[g] += lgi2[g];
		agi4[g] *= Cgxi[g];
	}
	// agj4
	for(size_t g=0; g<eg_size; ++g){
		agj4[g] = 0.0;
		for(size_t from_g=0; from_g<eg_size; ++from_g)
			agj4[g] += cal_bgk(g,from_g,Dgj,srgj,chigj,ssgj,vsfgj,keff)*agj2[from_g];
		agj4[g] += lgj2[g];
		agj4[g] *= Cgxj[g];
	}
	// agi1, agj1
	// block[0][0]
	{
	MAT *A = mat_create(2 * eg_size);
	for(size_t g=0; g<eg_size; ++g)
		for(size_t from_g=0; from_g<eg_size; ++from_g){
			double val = (-dgxi[g]) * (delta_func(g,from_g)+Fgxi[g]*cal_bgk(g,from_g,Dgi,srgi,chigi,ssgi,vsfgi,keff));
			mat_set(A, g, from_g, val);
		}
	// block[0][1]
	for(size_t g=0; g<eg_size; ++g)
		for(size_t from_g=0; from_g<eg_size; ++from_g){
			double val = dgxj[g] * (delta_func(g,from_g)+Fgxj[g]*cal_bgk(g,from_g,Dgj,srgj,chigj,ssgj,vsfgj,keff));
			mat_set(A, g, from_g+eg_size, val);
		}
	// block[1][0]
	for(size_t g=0; g<eg_size; ++g)
		for(size_t from_g=0; from_g<eg_size; ++from_g){
			double val = adfgi[g] * (delta_func(g,from_g)+Agxi[g]*cal_bgk(g,from_g,Dgi,srgi,chigi,ssgi,vsfgi,keff));
			mat_set(A, g+eg_size, from_g, val);
		}
	// block[1][1]
	for(size_t g=0; g<eg_size; ++g)
		for(size_t from_g=0; from_g<eg_size; ++from_g){
			double val = adfgj[g] * (delta_func(g,from_g)+Agxj[g]*cal_bgk(g,from_g,Dgj,srgj,chigj,ssgj,vsfgj,keff));
			mat_set(A, g+eg_size, from_g+eg_size, val);
		}
	VEC *b = vec_create(2 * eg_size);
	for(size_t g=0; g<eg_size; ++g){
		double val = dgxi[g] * (3.0*agi2[g] + Ggxi[g]*agi4[g] + Fgxi[g]*lgi1[g]) + dgxj[g] * (3.0*agj2[g] + Ggxj[g]*agj4[g] - Fgxj[g]*lgj1[g]);
		vec_set(b, g, val);
	}
	for(size_t g=0; g<eg_size; ++g){
		double val = (-adfgi[g]) * (phigi[g] + agi2[g] + agi4[g] + Agxi[g]*lgi1[g]) + adfgj[g] * (phigj[g] + agj2[g] + agj4[g] - Agxj[g]*lgj1[g]);
		vec_set(b, g+eg_size, val);
	}
	VEC *x = vec_create(2 * eg_size);
	LU_solve(x, A, b);
	for(size_t g=0; g<eg_size; ++g){
		agi1[g] = vec_get(x, g);
		agj1[g] = vec_get(x, g+eg_size);
	}
	vec_free(x);
	vec_free(b);
	mat_free(A);
	}
	// agi3
	for(size_t g=0; g<eg_size; ++g){
		agi3[g] = 0.0;
		for(size_t from_g=0; from_g<eg_size; ++from_g)
			agi3[g] += cal_bgk(g,from_g,Dgi,srgi,chigi,ssgi,vsfgi,keff)*agi1[from_g];
		agi3[g] += lgi1[g];
		agi3[g] *= Agxi[g];
	}
	// agj3
	for(size_t g=0; g<eg_size; ++g){
		agj3[g] = 0.0;
		for(size_t from_g=0; from_g<eg_size; ++from_g)
			agj3[g] += cal_bgk(g,from_g,Dgj,srgj,chigj,ssgj,vsfgj,keff)*agj1[from_g];
		agj3[g] += lgj1[g];
		agj3[g] *= Agxj[g];
	}
	for(size_t g=0; g<eg_size; ++g){
		J[g] = -dgxi[g]*(agi1[g] + 3.0*agi2[g] + Hgxi[g]*agi3[g] + Ggxi[g]*agi4[g]);
		#ifdef DEBUG
		double tmp = -dgxj[g]*(agj1[g] - 3.0*agj2[g] + Hgxj[g]*agj3[g] - Ggxj[g]*agj4[g]);
		if(fabs(J[g] - tmp) > 1e-6){
			fprintf(stderr, "SANM error.\n");
			exit(-1);
		}
		#endif
	}
	
	free(Agxi);
	free(Agxj);
	free(Cgxi);
	free(Cgxj);
	free(Egxi);
	free(Egxj);
	free(Fgxi);
	free(Fgxj);
	free(Ggxi);
	free(Ggxj);
	free(Hgxi);
	free(Hgxj);
	free(mgxi0);
	free(mgxj0);
	free(mgxi1);
	free(mgxj1);
	free(mgxi2);
	free(mgxj2);
	free(alpha_gxi);
	free(alpha_gxj);
	free(agi1);
	free(agj1);
	free(agi2);
	free(agj2);
	free(agi3);
	free(agj3);
	free(agi4);
	free(agj4);
	free(dgxi);
	free(dgxj);
}
