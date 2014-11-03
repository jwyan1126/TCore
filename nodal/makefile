NODAL = tnsol.o leak.o jcur.o
nodal.a($(NODAL)): $(NODAL)
	make -C SANM
	ar rc nodal.a $(NODAL)
tnsol.o: tnsol.h tnsol.c
	cc -c -g -Wall -std=c99 tnsol.h tnsol.c -DDEBUG
leak.o: leak.h leak.c
	cc -c -g -Wall -std=c99 leak.h leak.c -DDEBUG
jcur.o: jcur.h jcur.c
	cc -c -g -Wall -std=c99 jcur.h jcur.c -DDEBUG
clean:
	make clean -C SANM
	rm *.o
	rm *.gch
	rm *.a
