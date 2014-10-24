#ifndef MTRLLIB_H
#define MTRLLIB_H

#include"mtrl.h"
#include"list.h"
#include<stdio.h>

typedef LIST MTRLLIB;

MTRLLIB *mtrllib_create();

void mtrllib_free(MTRLLIB *mlib);

size_t mtrllib_get_size(const MTRLLIB *mlib);

void mtrllib_add(MTRLLIB *mlib, MTRL *mtrl);

void mtrllib_remove_fromid(MTRLLIB *mlib, int mtrl_id);

MTRL *mtrllib_get_fromid(const MTRLLIB *mlib, int mtrl_id);

void mtrllib_fprintf(const MTRLLIB *mlib, FILE *stream);

#endif
