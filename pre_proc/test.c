#include<stdio.h>
#include"mtrl.h"
#include"mtrllib.h"
#include"input.h"
#include"sconf.h"
#include<stdlib.h>

int main()
{
	INPUT *input = input_create(NULL);
	SCONF *sconf = sconf_create(input);
	sconf_fprintf(sconf, stdout);
	sconf_free(sconf);
	input_free(input);
	return 0;
}
