#include"sanm.h"
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"../../algebra/mat.h"
#include"../../algebra/vec.h"
#include"../../algebra/ksp.h"
#include"cal_bgk.h"
#include<limits.h>

void sanm_left(TNSOL *tn)
{
	size_t eg_size = tn->eg_size;
	int bdy = tn->bdy;
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
	double *agj4 = malloc(eg_size * sizeof(double));
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
	// agj4
	for(size_t g=0; g<eg_size; ++g){
		agj4[g] = 0.0;
		for(size_t from_g=0; from_g<eg_size; ++from_g)
			agj4[g] += cal_bgk(g,from_g,Dgj,srgj,chigj,ssgj,vsfgj,keff)*agj2[from_g];
		agj4[g] += lgj2[g];
		agj4[g] *= Cgxj[g];
	}
	// agj1
	if(bdy == 0){
		for(size_t g=0; g<eg_size; ++g)
			J[g] = 0.0;
		return;
	}
	else if(bdy == 1){
		MAT *A = mat_create(eg_size);
		for(size_t g=0; g<eg_size; ++g)
			for(size_t from_g=0; from_g<eg_size; ++from_g){
				double val = delta_func(g,from_g) + Agxj[g]*cal_bgk(g,from_g,Dgj,srgj,chigj,ssgj,vsfgj,keff);
				mat_set(A, g, from_g, val);
			}
		VEC *b = vec_create(eg_size);
		for(size_t g=0; g<eg_size; ++g){
			double val = phigj[g] + agj2[g] + agj4[g] - Agxj[g]*lgj1[g];
			vec_set(b, g, val);
		}
		VEC *x = vec_create(eg_size);
		LU_solve(x, A, b);
		for(size_t g=0; g<eg_size; ++g)
			agj1[g] = vec_get(x, g);
		vec_free(x);
		vec_free(b);
		mat_free(A);
	}
	else if(bdy == 2){
		MAT *A = mat_create(eg_size);
		for(size_t g=0; g<eg_size; ++g)
			for(size_t from_g=0; from_g<eg_size; ++from_g){
				double val = (2.0*dgxj[g]+1)*delta_func(g,from_g) + (Agxj[g]+2.0*dgxj[g]*Fgxj[g])*cal_bgk(g,from_g,Dgj,srgj,chigj,ssgj,vsfgj,keff);
				mat_set(A, g, from_g, val);
			}
		VEC *b = vec_create(eg_size);
		for(size_t g=0; g<eg_size; ++g){
			double val = phigj[g] + (1+6.0*dgxj[g])*agj2[g] + (1+2.0*dgxj[g]*Ggxj[g])*agj4[g] - (Agxj[g]+2.0*dgxj[g]*Fgxj[g])*lgj1[g];
			vec_set(b, g, val);
		}
		VEC *x = vec_create(eg_size);
		LU_solve(x, A, b);
		for(size_t g=0; g<eg_size; ++g)
			agj1[g] = vec_get(x, g);
		vec_free(x);
		vec_free(b);
		mat_free(A);
	}
	else{ fprintf(stderr, "SANM error.\n"); exit(-1); }
	// agj3
	for(size_t g=0; g<eg_size; ++g){
		agj3[g] = 0.0;
		for(size_t from_g=0; from_g<eg_size; ++from_g)
			agj3[g] += cal_bgk(g,from_g,Dgj,srgj,chigj,ssgj,vsfgj,keff)*agj1[from_g];
		agj3[g] += lgj1[g];
		agj3[g] *= Agxj[g];
	}
	for(size_t g=0; g<eg_size; ++g)
		J[g] = -dgxj[g]*(agj1[g] - 3.0*agj2[g] + Hgxj[g]*agj3[g] - Ggxj[g]*agj4[g]);

	free(Agxj);
	free(Cgxj);
	free(Egxj);
	free(Fgxj);
	free(Ggxj);
	free(Hgxj);
	free(mgxj0);
	free(mgxj1);
	free(mgxj2);
	free(alpha_gxj);
	free(agj1);
	free(agj2);
	free(agj3);
	free(agj4);
	free(dgxj);
}
