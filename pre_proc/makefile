test: test.c list.o mtrl.o mtrllib.o
	cc -g -Wall -std=c99 test.c list.o mtrl.o mtrllib.o
list.o: list.h list.c
	cc -g -c -Wall -std=c99 list.h list.c -DDEBUG
mtrl.o: mtrl.h mtrl.c
	cc -g -c -Wall -std=c99 mtrl.h mtrl.c -DDEBUG
mtrllib.o: mtrllib.h mtrllib.c
	cc -g -c -Wall -std=c99 mtrllib.h mtrllib.c -DDEBUG
clean:
	rm *.o
	rm *.gch
