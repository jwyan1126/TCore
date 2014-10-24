#include<stdlib.h>
#include<stdio.h>
#include"rect_mapper.h"
#include"input.h"
#include"sconf.h"
#include"msh.h"

//size_t EG_SIZE; // num of energy groups
//size_t XM_SIZE; // max num of meshes in x direction
//size_t YM_SIZE; // max num of meshes in y direction
//size_t ZM_SIZE; // max num of meshes in z direction
//size_t RT_SIZE; // num of meshes in xyz direction(3D)
//RECT_MAPPER *MAPPER;

int main(int argc, char *argv[])
{
	INPUT *input = input_create(NULL);
	SCONF *steady_conf = sconf_create(input);
	MAPPER = rect_mapper_create();
	MSH *msh = msh_create(steady_conf);
	
	msh_free(msh);
	rect_mapper_free(MAPPER);
	sconf_free(steady_conf);
	input_free(input);
	return EXIT_SUCCESS;
}
