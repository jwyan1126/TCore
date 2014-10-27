#ifndef STEADY_SOLVER_H
#define STEADY_SOLVER_H

#include"pre_proc/sconf.h"
#include"pre_proc/edat.h"
#include"algebra/mat.h"

void steady_solver(SSOL *ssol, const SCONF *sconf, const MESH *mesh);

void cal_DFDM(EDAT *DFDM, const *sconf, const *mesh);

void cal_M(MAT *M, const EDAT *DFDM, const EDAT *DNOD);

void cal_S(MAT *S, const *sconf, const *mesh);

void cal_F(MAT *F, const (sconf, const *mesh);

#endif
