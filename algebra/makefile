algebra.a(vec.o mat.o mvop.o bicgstab.o gauss_seidel.o eigen.o): vec.o mat.o mvop.o bicgstab.o gauss_seidel.o eigen.o
	ar rc algebra.a vec.o mat.o mvop.o bicgstab.o gauss_seidel.o eigen.o
vec.o: vec.h vec.c
	cc -g -c -Wall -std=c99 vec.h vec.c -DDEBUG
mat.o: mat.h mat.c
	cc -g -c -Wall -std=c99 mat.h mat.c -DDEBUG
mvop.o: mvop.h mvop.c
	cc -g -c -Wall -std=c99 mvop.h mvop.c -fopenmp -DDEBUG
bicgstab.o: ksp.h bicgstab.c
	cc -g -c -Wall -std=c99 ksp.h bicgstab.c -DDEBUG
gauss_seidel.o: ksp.h gauss_seidel.c
	cc -g -c -Wall -std=c99 ksp.h gauss_seidel.c -DDEBUG
eigen.o: eigen.h eigen.c
	cc -g -c -Wall -std=c99 eigen.h eigen.c -DDEBUG
clean:
	rm *.o
	rm *.a
	rm *.gch
