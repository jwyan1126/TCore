FILES = main.c steady_solver.o transient_solver.o ssol.o tsol.o cal_DFDM.o cal_DNOD.o flux.o pcs.o control_rod.o cal_s.o cal_f.o cal_m.o power.o adjust_vsf.o
main: $(FILES)
	make -C algebra
	make -C pre_proc
	make -C nodal
	cc -g -Wall -std=c99 $(FILES) pre_proc/pre_proc.a algebra/algebra.a nodal/nodal.a nodal/SANM/SANM.a -fopenmp -lm -DDEBUG
steady_solver.o: steady_solver.h steady_solver.c
	cc -g -c -Wall -std=c99 steady_solver.h steady_solver.c -DDEBUG
transient_solver.o: transient_solver.h transient_solver.c
	cc -g -c -Wall -std=c99 transient_solver.h transient_solver.c -DDEBUG
cal_s.o: steady_solver.h cal_s.c
	cc -g -c -Wall -std=c99 steady_solver.h cal_s.c -DDEBUG
cal_m.o: steady_solver.h cal_m.c
	cc -g -c -Wall -std=c99 steady_solver.h cal_m.c -DDEBUG
cal_f.o: steady_solver.h cal_f.c
	cc -g -c -Wall -std=c99 steady_solver.h cal_f.c -DDEBUG
adjust_vsf.o: steady_solver.h adjust_vsf.c
	cc -g -c -Wall -std=c99 steady_solver.h adjust_vsf.c -DDEBUG
cal_DFDM.o: cal_DFDM.c steady_solver.h
	cc -g -c -Wall -std=c99 cal_DFDM.c steady_solver.h -DDEBUG
cal_DNOD.o: cal_DNOD.c steady_solver.h
	cc -g -c -Wall -std=c99 cal_DNOD.c steady_solver.h -DDEBUG
ssol.o: ssol.h ssol.c
	cc -g -c -Wall -std=c99 ssol.h ssol.c -DDEBUG
tsol.o: tsol.h tsol.c
	cc -g -c -Wall -std=c99 tsol.h tsol.c -DDEBUG
flux.o: flux.h flux.c
	cc -g -c -Wall -std=c99 flux.h flux.c -DDEBUG
pcs.o: pcs.h pcs.c
	cc -g -c -Wall -std=c99 pcs.h pcs.c -DDEBUG
control_rod.o: control_rod.h control_rod.c
	cc -g -c -Wall -std=c99 control_rod.h control_rod.c -DDEBUG
power.o: power.h power.c
	cc -g -c -Wall -std=c99 power.h power.c -DDEBUG
clean:
	make clean -C algebra
	make clean -C pre_proc
	make clean -C nodal
	rm *.gch
	rm *.o
	rm a.out
