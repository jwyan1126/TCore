#include<stdio.h>
#include"mtrl.h"
#include"mtrllib.h"
#include"input.h"
#include"mapper.h"
#include"sconf.h"
#include<stdlib.h>
#include"mesh.h"
#include"cdat.h"
#include"edat.h"

int main()
{
	INPUT *input = input_create(NULL);
	SCONF *sconf = sconf_create(input);
	MAPPER *mapper = mapper_create(sconf);
	MESH *mesh = mesh_create(sconf);
	EDAT4 *edat = edat4_create(sconf);
	edat4_set_rand(edat);
	edat4_fprintf(edat, 0, 5, 5, 0, stdout);
	
	edat4_free(edat);
	mesh_free(mesh);
	mapper_free(mapper);
	sconf_free(sconf);
	input_free(input);
	return 0;
}
