#include<stdio.h>
#include"pre_proc/mtrl.h"
#include"pre_proc/mtrllib.h"
#include"pre_proc/input.h"
#include"pre_proc/mapper.h"
#include"pre_proc/sconf.h"
#include<stdlib.h>
#include"pre_proc/mesh.h"
#include"ssol.h"
#include"steady_solver.h"
#include"algebra/krylov.h"

MATSOLVER mat_solver;

int main()
{
	mat_solver = bicgstab;
	INPUT *input = input_create(NULL);
	SCONF *sconf = sconf_create(input);
	MAPPER *mapper = mapper_create(sconf);
	MESH *mesh = mesh_create(sconf);
	SSOL *ssol = ssol_create(sconf->eg_size, sconf->rt_mesh_size);
	steady_solver(ssol, sconf, mapper, mesh);
	printf("%g\n", ssol->keff);
	
	ssol_free(ssol);
	mesh_free(mesh);
	mapper_free(mapper);
	sconf_free(sconf);
	input_free(input);
	return 0;
}
