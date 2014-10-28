main: main.c steady_solver.o ssol.o cal_DFDM.o
	make -C algebra
	make -C pre_proc
	cc -g -Wall -std=c99 main.c ssol.o steady_solver.o cal_DFDM.o pre_proc/pre_proc.a algebra/algebra.a -lm -DDEBUG
steady_solver.o: steady_solver.h steady_solver.c
	cc -g -c -Wall -std=c99 steady_solver.h steady_solver.c -DDEBUG
cal_DFDM.o: cal_DFDM.c steady_solver.h
	cc -g -c -Wall -std=c99 cal_DFDM.c steady_solver.h -DDEBUG
ssol.o: ssol.h ssol.c
	cc -g -c -Wall -std=c99 ssol.h ssol.c -DDEBUG
clean:
	make clean -C algebra
	make clean -C pre_proc
	rm a.out
	rm *.gch
	rm *.o
