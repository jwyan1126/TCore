#ifndef JCUR_H
#define JCUR_H

#include"../pre_proc/mesh.h"
#include"../pre_proc/edat.h"
#include"leak.h"
#include"../ssol.h"

void cal_jcur(EDAT4 *jcur, const MESH *mesh, const LEAK *leak, const SSOL *ssol);

#endif
