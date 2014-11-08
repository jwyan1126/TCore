#include<stdio.h>
#include"pre_proc/mtrl.h"
#include"pre_proc/mtrllib.h"
#include"pre_proc/input.h"
#include"pre_proc/mapper.h"
#include"pre_proc/sconf.h"
#include"pre_proc/tconf.h"
#include<stdlib.h>
#include"pre_proc/mesh.h"
#include"ssol.h"
#include"tsol.h"
#include"steady_solver.h"
#include"transient_solver.h"
#include"algebra/ksp.h"
#include"control_rod.h"
#include"power.h"

MATSOLVER mat_solver;

int main()
{
	mat_solver = bicgstab;
	INPUT *input = input_create(NULL);
	SCONF *sconf = sconf_create(input);
	TCONF *tconf = tconf_create(input);
	MAPPER *mapper = mapper_create(sconf);
	MESH *mesh = mesh_create(sconf, mapper);
	FLUX *flux = flux_create(mapper);
	PCS *pcs = pcs_create(tconf->pcs_size, mapper);
	double keff = steady_state_cal(flux, sconf, mapper, mesh);
	flux_normalize(flux);
	printf("Steady cal. finished. keff = %g\n", keff);
	adjust_vsf(mesh->vsf, keff);
	pcs_init(pcs, flux, mesh, tconf);
	int steps = tconf->steps;
	double tau = tconf->tau;
	CDAT4 *vsf_lst = cdat4_create(mapper);
	double cur_time = 0.0;
	for(int step = 0; step < steps; ++step){
		cur_time += tau;
		cdat4_copy(vsf_lst, mesh->vsf);
		cross_section_update(mesh);
		next_step_cal(flux, pcs, tconf, mapper, mesh, vsf_lst);
		printf("%g\n",flux_sumup(flux));
	}

	pcs_free(pcs);
	flux_free(flux);
	cdat4_free(vsf_lst);
	tconf_free(tconf);
	mesh_free(mesh);
	mapper_free(mapper);
	sconf_free(sconf);
	input_free(input);
	return 0;
}
