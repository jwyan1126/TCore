ALGEBRA = vec.o mat.o mvop.o ksp.o eigen.o
algebra.a($(ALGEBRA)): $(ALGEBRA) 
	ar rc algebra.a $(ALGEBRA)
vec.o: vec.h vec.c
	cc -g -c -Wall -std=c99 vec.h vec.c -DDEBUG
mat.o: mat.h mat.c
	cc -g -c -Wall -std=c99 mat.h mat.c -DDEBUG
mvop.o: mvop.h mvop.c
	cc -g -c -Wall -std=c99 mvop.h mvop.c -fopenmp -DDEBUG
ksp.o: ksp.h ksp.c
	cc -g -c -Wall -std=c99 ksp.h ksp.c -DDEBUG
eigen.o: eigen.h eigen.c
	cc -g -c -Wall -std=c99 eigen.h eigen.c -DDEBUG
clean:
	rm *.o
	rm *.a 
	rm *.gch
