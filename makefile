FILES = main.c steady_solver.o ssol.o cal_DFDM.o cal_DNOD.o flux.o
main: $(FILES)
	make -C algebra
	make -C pre_proc
	make -C nodal
	cc -g -Wall -std=c99 $(FILES) pre_proc/pre_proc.a algebra/algebra.a nodal/nodal.a nodal/SANM/SANM.a -fopenmp -lm -DDEBUG
steady_solver.o: steady_solver.h steady_solver.c
	cc -g -c -Wall -std=c99 steady_solver.h steady_solver.c -DDEBUG
cal_DFDM.o: cal_DFDM.c steady_solver.h
	cc -g -c -Wall -std=c99 cal_DFDM.c steady_solver.h -DDEBUG
cal_DNOD.o: cal_DNOD.c steady_solver.h
	cc -g -c -Wall -std=c99 cal_DNOD.c steady_solver.h -DDEBUG
ssol.o: ssol.h ssol.c
	cc -g -c -Wall -std=c99 ssol.h ssol.c -DDEBUG
flux.o: flux.h flux.c
	cc -g -c -Wall -std=c99 flux.h flux.c -DDEBUG
clean:
	make clean -C algebra
	make clean -C pre_proc
	make clean -C nodal
	rm *.gch
	rm *.o
	rm a.out
