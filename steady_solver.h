#ifndef STEADY_SOLVER_H
#define STEADY_SOLVER_H

#include"pre_proc/sconf.h"
#include"pre_proc/edat.h"
#include"algebra/mat.h"
#include"ssol.h"
#include"pre_proc/mesh.h"

void steady_solver(SSOL *ssol, SCONF *sconf, MAPPER *mapper, const MESH *mesh);

void cal_DFDM(EDAT4 *DFDM, const MESH *mesh);

void cal_DNOD(EDAT4 *DNOD, const MESH *mesh, const EDAT4 *DFDM, const EDAT4 *Jn, const SSOL *ssol);

void cal_M(MAT *M, EDAT4 *DFDM, EDAT4 *DNOD, const MAPPER *mapper, const MESH *mesh);

void cal_S(MAT *S, const SCONF *sconf, const MAPPER *mapper, const MESH *mesh);

void cal_F(MAT *F, const SCONF *sconf, const MAPPER *mapper, const MESH *mesh);

#endif
