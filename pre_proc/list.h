#ifndef LIST_H
#define LIST_H

#include<stddef.h>
#include<stdio.h>
#include"mtrl.h"

typedef MTRL *LISTTYPE;

struct NODE
{
	LISTTYPE data;
	struct NODE *next;
};

typedef struct
{
	size_t size;
	struct NODE *head;
} LIST;

LIST *list_create();

void list_free(LIST *list);

size_t list_len(const LIST *list);

LISTTYPE list_get_val(const LIST *list, size_t index);

void list_fprintf(const LIST *list, FILE *stream, char *to_string(LISTTYPE listtype));

// Insert element at the beginning
void list_push(LIST *list, LISTTYPE val);

// Remove element at the beginning
LISTTYPE list_pop(LIST *list);

void list_insert(LIST *list, size_t index, LISTTYPE val);

LISTTYPE list_remove(LIST *list, size_t index);

void list_clear(LIST *list);

#endif
