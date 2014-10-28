algebra.a(vec.o mat.o mvop.o bicgstab.o eigen.o): vec.o mat.o mvop.o bicgstab.o eigen.o
	ar rc algebra.a vec.o mat.o mvop.o bicgstab.o eigen.o
vec.o: vec.h vec.c
	cc -g -c -Wall -std=c99 vec.h vec.c -DDEBUG
mat.o: mat.h mat.c
	cc -g -c -Wall -std=c99 mat.h mat.c -DDEBUG
mvop.o: mvop.h mvop.c
	cc -g -c -Wall -std=c99 mvop.h mvop.c -DDEBUG
bicgstab.o: krylov.h bicgstab.c
	cc -g -c -Wall -std=c99 krylov.h bicgstab.c -DDEBUG
eigen.o: eigen.h eigen.c
	cc -g -c -Wall -std=c99 eigen.h eigen.c -DDEBUG
clean:
	rm *.o
	rm *.a
