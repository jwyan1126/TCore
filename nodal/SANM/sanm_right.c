#include"sanm.h"
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"../../algebra/mat.h"
#include"../../algebra/vec.h"
#include"../../algebra/ksp.h"
#include"cal_bgk.h"
#include<limits.h>

void sanm_right(TNSOL *tn)
{
	size_t eg_size = tn->eg_size;
	int bdy = tn->bdy;
	double dui = tn->dui;
	double *Dgi = tn->Dgi;
	double *vsfgi = tn->vsfgi;
	double *phigi = tn->phigi;
	double *srgi = tn->srgi;
	double **ssgi = tn->ssgi;
	double *chigi = tn->chigi;
	double *lgi0 = tn->lgi0;
	double *lgi1 = tn->lgi1;
	double *lgi2 = tn->lgi2;
	double keff = tn->keff;
	double *J = tn->J;
	double *dgxi = malloc(eg_size * sizeof(double));
	double *agi1 = malloc(eg_size * sizeof(double));
	double *agi2 = malloc(eg_size * sizeof(double));
	double *agi3 = malloc(eg_size * sizeof(double));
	double *agi4 = malloc(eg_size * sizeof(double));
	double *alpha_gxi = malloc(eg_size * sizeof(double));
	double *mgxi0 = malloc(eg_size * sizeof(double));
	double *mgxi1 = malloc(eg_size * sizeof(double));
	double *mgxi2 = malloc(eg_size * sizeof(double));
	double *Agxi = malloc(eg_size * sizeof(double));
	double *Cgxi = malloc(eg_size * sizeof(double));
	double *Egxi = malloc(eg_size * sizeof(double));
	double *Fgxi = malloc(eg_size * sizeof(double));
	double *Ggxi = malloc(eg_size * sizeof(double));
	double *Hgxi = malloc(eg_size * sizeof(double));

	for(size_t g=0; g<eg_size; ++g){
		dgxi[g] = 2.0*Dgi[g] / dui;
		alpha_gxi[g] = sqrt(srgi[g]/Dgi[g])*dui/2.0;
    		mgxi0[g] = sinh(alpha_gxi[g])/alpha_gxi[g];
    		mgxi1[g] = 3.0*(alpha_gxi[g]*cosh(alpha_gxi[g])-sinh(alpha_gxi[g]))/(alpha_gxi[g]*alpha_gxi[g]); 
    		mgxi2[g] = 5.0*(alpha_gxi[g]*alpha_gxi[g]*sinh(alpha_gxi[g])-3.0*alpha_gxi[g]*cosh(alpha_gxi[g])+3.0*sinh(alpha_gxi[g]))/(alpha_gxi[g]*alpha_gxi[g]*alpha_gxi[g]);
    		Agxi[g] = (sinh(alpha_gxi[g])-mgxi1[g])/(srgi[g]*mgxi1[g]);
    		Cgxi[g] = (cosh(alpha_gxi[g])-mgxi0[g]-mgxi2[g])/(srgi[g]*mgxi2[g]);
    		Egxi[g] = (mgxi0[g]/mgxi2[g])-3.0/(alpha_gxi[g]*alpha_gxi[g]);
    		Fgxi[g] = (alpha_gxi[g]*cosh(alpha_gxi[g])-mgxi1[g])/(srgi[g]*mgxi1[g]);
    		Ggxi[g] = (alpha_gxi[g]*sinh(alpha_gxi[g])-3.0*mgxi2[g])/(cosh(alpha_gxi[g])-mgxi0[g]-mgxi2[g]);
    		Hgxi[g] = (alpha_gxi[g]*cosh(alpha_gxi[g])-mgxi1[g])/(sinh(alpha_gxi[g])-mgxi1[g]);
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
	LU_solve(x, A, b);
	for(size_t g=0; g<eg_size; ++g)
		agi2[g] = vec_get(x, g);
	vec_free(x);
	vec_free(b);
	mat_free(A);
	// agi4
	for(size_t g=0; g<eg_size; ++g){
		agi4[g] = 0.0;
		for(size_t from_g=0; from_g<eg_size; ++from_g)
			agi4[g] += cal_bgk(g,from_g,Dgi,srgi,chigi,ssgi,vsfgi,keff)*agi2[from_g];
		agi4[g] += lgi2[g];
		agi4[g] *= Cgxi[g];
	}
	// aji1
	if(bdy == 0){
		for(size_t g=0; g<eg_size; ++g)
			J[g] = 0.0;
		return;
	}
	else if(bdy == 1){
		MAT *A = mat_create(eg_size);
		for(size_t g=0; g<eg_size; ++g)
			for(size_t from_g=0; from_g<eg_size; ++from_g){
				double val = delta_func(g,from_g) + Agxi[g]*cal_bgk(g,from_g,Dgi,srgi,chigi,ssgi,vsfgi,keff);
				mat_set(A, g, from_g, val);
			}
		VEC *b = vec_create(eg_size);
		for(size_t g=0; g<eg_size; ++g){
			double val = -phigi[g] - agi2[g] - agi4[g] - Agxi[g]*lgi1[g];
			vec_set(b, g, val);
		}
		VEC *x = vec_create(eg_size);
		LU_solve(x, A, b);
		for(size_t g=0; g<eg_size; ++g)
			agi1[g] = vec_get(x, g);
		vec_free(x);
		vec_free(b);
		mat_free(A);
	}
	else if(bdy == 2){
		MAT *A = mat_create(eg_size);
		for(size_t g=0; g<eg_size; ++g)
			for(size_t from_g=0; from_g<eg_size; ++from_g){
				double val = (2.0*dgxi[g]+1)*delta_func(g,from_g) + (Agxi[g]+2.0*dgxi[g]*Fgxi[g])*cal_bgk(g,from_g,Dgi,srgi,chigi,ssgi,vsfgi,keff);
				mat_set(A, g, from_g, val);
			}
		VEC *b = vec_create(eg_size);
		for(size_t g=0; g<eg_size; ++g){
			double val = -phigi[g] - (1+6.0*dgxi[g])*agi2[g] - (1+2.0*dgxi[g]*Ggxi[g])*agi4[g] - (Agxi[g]+2.0*dgxi[g]*Fgxi[g])*lgi1[g];
			vec_set(b, g, val);
		}
		VEC *x = vec_create(eg_size);
		LU_solve(x, A, b);
		for(size_t g=0; g<eg_size; ++g)
			agi1[g] = vec_get(x, g);
		vec_free(x);
		vec_free(b);
		mat_free(A);
	}
	else{ fprintf(stderr, "SANM error.\n"); exit(-1); }
	// agi3
	for(size_t g=0; g<eg_size; ++g){
		agi3[g] = 0.0;
		for(size_t from_g=0; from_g<eg_size; ++from_g)
			agi3[g] += cal_bgk(g,from_g,Dgi,srgi,chigi,ssgi,vsfgi,keff)*agi1[from_g];
		agi3[g] += lgi1[g];
		agi3[g] *= Agxi[g];
	}
	for(size_t g=0; g<eg_size; ++g)
		J[g] = -dgxi[g]*(agi1[g] + 3.0*agi2[g] + Hgxi[g]*agi3[g] + Ggxi[g]*agi4[g]);

	free(Agxi);
	free(Cgxi);
	free(Egxi);
	free(Fgxi);
	free(Ggxi);
	free(Hgxi);
	free(mgxi0);
	free(mgxi1);
	free(mgxi2);
	free(alpha_gxi);
	free(agi1);
	free(agi2);
	free(agi3);
	free(agi4);
	free(dgxi);
}
