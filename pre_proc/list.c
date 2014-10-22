#include"list.h"
#include<stdlib.h>

LIST *list_create()
{
	LIST *list = malloc(sizeof(LIST));
	list->size = 0;
	list->head = NULL;
	return list;
}

void list_free(LIST *list)
{
	struct NODE *cur_node = list->head;
	struct NODE *next_node;
	while(cur_node != NULL){
		next_node = cur_node->next;
		free(cur_node);
		cur_node = next_node;
	}
	free(list);
}

inline size_t list_len(const LIST *list)
{
	return list->size;
}

struct NODE *list_get_node(const LIST *list, size_t index)
{
	#ifdef DEBUG
	if(index >= list_len(list)){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	struct NODE *node = list->head;
	for(size_t i=0; i<index; ++i)
		node = node->next;
	return node;
}

LISTTYPE list_get_val(const LIST *list, size_t index)
{
	#ifdef DEBUG
	if(index >= list_len(list)){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	struct NODE *node = list->head;
	for(size_t i=0; i<index; ++i)
		node = node->next;
	return node->data;
}

void list_fprintf(const LIST *list, FILE *stream, char *to_string(LISTTYPE listtype))
{
	struct NODE *cur_node = list->head;
	while(cur_node!=NULL){
		fprintf(stream, "%s\n", to_string(cur_node->data));
		cur_node = cur_node->next;
	}
	fprintf(stream, "\n");
}

// Insert element at the beginning
void list_push(LIST *list, LISTTYPE val)
{
	struct NODE *node = malloc(sizeof(struct NODE));
	node->data = val;
	node->next = list->head;
	list->head = node;
	list->size++;
}

// Remove element at the beginning
LISTTYPE list_pop(LIST *list)
{
	#ifdef DEBUG
	if(!list_len(list)){
		fprintf(stderr, "List have no element.\n");
		exit(-1);
	}
	#endif
	struct NODE *node = list->head->next;
	LISTTYPE val = list->head->data;
	free(list->head);
	list->head = node;
	list->size--;
	return val;
}

void list_insert(LIST *list, size_t index, LISTTYPE val)
{
	#ifdef DEBUG
	if(index > list_len(list)){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	struct NODE *node = malloc(sizeof(struct NODE));
	node->data = val;
	if(index == 0)
		list_push(list, val);
	else{
		struct NODE *former_node = list_get_node(list, index-1);
		node->next = former_node->next;
		former_node->next = node;
	}
	list->size++;
}

LISTTYPE list_remove(LIST *list, size_t index)
{
	#ifdef DEBUG
	if(index >= list_len(list)){
		fprintf(stderr, "Index out of range.\n");
		exit(-1);
	}
	#endif
	LISTTYPE val;
	if(index == 0)
		val = list_pop(list);
	else{
		struct NODE *former_node = list_get_node(list, index-1);
		struct NODE *node = former_node->next;
		former_node->next = node->next;
		val = node->data;
		free(node);
	}
	list->size--;
	return val;
}

void list_clear(LIST *list)
{
	struct NODE *cur_node = list->head;
	struct NODE *next_node;
	while(cur_node != NULL){
		next_node = cur_node->next;
		free(cur_node);
		cur_node = next_node;
		list->size--;
	}
	list->head = NULL;
}
