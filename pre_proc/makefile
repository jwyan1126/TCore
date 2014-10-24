test: test.c list.o mtrl.o mtrllib.o input.o sconf.o mapper.o msh.o
	cc -g -Wall -std=c99 test.c list.o mtrl.o mtrllib.o input.o sconf.o mapper.o msh.o -DDEBUG
input.o: input.h input.c
	cc -g -c -Wall -std=c99 input.h input.c -DDEBUG
sconf.o: sconf.h sconf.c
	cc -g -c -Wall -std=c99 sconf.h sconf.c -DDEBUG
mapper.o: mapper.h mapper.c
	cc -g -c -Wall -std=c99 mapper.h mapper.c -DDEBUG
list.o: list.h list.c
	cc -g -c -Wall -std=c99 list.h list.c -DDEBUG
mtrl.o: mtrl.h mtrl.c
	cc -g -c -Wall -std=c99 mtrl.h mtrl.c -DDEBUG
mtrllib.o: mtrllib.h mtrllib.c
	cc -g -c -Wall -std=c99 mtrllib.h mtrllib.c -DDEBUG
msh.o: msh.h msh.c
	cc -g -c -Wall -std=c99 msh.h msh.c -DDEBUG
clean:
	rm *.o
	rm *.gch
	rm a.out
