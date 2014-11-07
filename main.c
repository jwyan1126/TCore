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
	MAPPER *mapper = mapper_create(sconf);
	MESH *mesh = mesh_create(sconf, mapper);
	SSOL *ssol = ssol_create(mapper);
	POWER *power = power_create(mapper);
	steady_state_cal(ssol, sconf, mapper, mesh);
	printf("steady calculation finished.\n");
	printf("keff = %g\n", ssol->keff);
	adjust_vsf(mesh->vsf, ssol->keff);
	flux_normalize(ssol->flux);
	TCONF *tconf = tconf_create(input);
	TSOL *tsol = tsol_create(tconf, mesh, ssol);
	int steps = tconf->steps;
	for(int step = 0; step < steps; ++step){
		cross_section_update(mesh);
		next_step_cal(tsol, tconf, mapper, mesh);
		//power_cal(power, tsol->flux, mesh->vsf);
		printf("%g\n", flux_sumup(tsol->flux));
	}
	
	power_free(power);
	tsol_free(tsol);
	tconf_free(tconf);
	ssol_free(ssol);
	mesh_free(mesh);
	mapper_free(mapper);
	sconf_free(sconf);
	input_free(input);
	return 0;
}
