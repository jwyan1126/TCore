#ifndef TRANSIENT_SOLVER_H
#define TRANSIENT_SOLVER_H

#include"algebra/vec.h"
#include"tsol.h"
#include"pre_proc/mapper.h"
#include"pre_proc/mesh.h"
#include"pre_proc/tconf.h"
#include"flux.h"
#include"pre_proc/cdat.h"
#include"pcs.h"
#include"steady_solver.h"
#include"algebra/ksp.h"

extern MATSOLVER mat_solver;

void next_step_cal(TSOL *tsol, TCONF *tconf, MAPPER *mapper, MESH *mesh);

void src_rvs_cal(VEC *src, TCONF *tconf, MESH *mesh, CDAT4 *vsf_lst, FLUX *flux, FLUX *flux_lst, PCS *pcs, PCS *pcs_lst);

void sr_rvs_cal(CDAT4 *sr_rvs, TCONF *tconf, MESH *mesh);

void pcs_cal(PCS *pcs, FLUX *flux, FLUX *flux_lst, CDAT4 *vsf, CDAT4 *vsf_lst, TCONF *tconf);

#endif
