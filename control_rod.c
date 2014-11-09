#include"control_rod.h"

void cross_section_update(MESH *mesh, double cur_time, SCONF *sconf)
{
	MBLOCK mblock = sconf_get_mblock(sconf, 0,0,0);
	size_t start_x = mblock.start_x;
	size_t end_x = mblock.end_x;
	size_t start_y = mblock.start_y;
	size_t end_y = mblock.end_y;
	size_t start_z = mblock.start_z;
	size_t end_z = mblock.end_z;
	for(size_t i=start_x; i<=end_x; ++i)
		for(size_t j=start_y; j<=end_y; ++j)
			for(size_t k=start_z; k<=end_z; ++k){
				mesh->sr->data[k][j][i][1] = 0.18*(1.0-0.01*cur_time);
			}
}
