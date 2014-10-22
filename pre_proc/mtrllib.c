#include"mtrllib.h"
#include<stdlib.h>

MTRLLIB *mtrllib_create()
{
	return list_create();
}

void mtrllib_free(MTRLLIB *mlib)
{
	list_free(mlib);
}

size_t mtrllib_get_size(const MTRLLIB *mlib)
{
	return mlib->size;
}

void mtrllib_add(MTRLLIB *mlib, MTRL *mtrl)
{
	list_push(mlib, mtrl);
}

void mtrllib_remove_fromid(MTRLLIB *mlib, int mtrl_id)
{
	size_t index = 0;
	struct NODE *cur_node = mlib->head;
	while(cur_node != NULL){
		if(cur_node->data->mtrl_id == mtrl_id)
			break;
		cur_node = cur_node->next;
		index++;
	}
	if(index < mlib->size)
		list_remove(mlib, index);
	else{
		fprintf(stderr, "mtrl_id missing.\n");
		exit(-1);
	}
}

MTRL *mtrllib_get_fromid(const MTRLLIB *mlib, int mtrl_id)
{
	struct NODE *cur_node = mlib->head;
	while(cur_node != NULL){
		if(cur_node->data->mtrl_id == mtrl_id)
			return cur_node->data;
		cur_node = cur_node->next;
	}
	return NULL;
}
