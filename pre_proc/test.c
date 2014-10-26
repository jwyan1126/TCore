#include<stdio.h>
#include"mtrl.h"
#include"mtrllib.h"
#include"input.h"
#include"mapper.h"
#include"sconf.h"
#include<stdlib.h>
#include"mesh.h"

int main()
{
	INPUT *input = input_create(NULL);
	SCONF *sconf = sconf_create(input);
	MAPPER *mapper = mapper_create(sconf);
	MESH *mesh = mesh_create(sconf);

	mesh_free(mesh);
	mapper_free(mapper);
	sconf_free(sconf);
	input_free(input);
	return 0;
}
