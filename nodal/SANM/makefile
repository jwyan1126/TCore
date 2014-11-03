SANM = sanm_inner.o sanm_left.o sanm_right.o cal_bgk.o
SANM.a($(SANM)): $(SANM)
	ar rc SANM.a $(SANM)
sanm_inner.o: sanm.h sanm_inner.c
	cc -c -g -Wall -std=c99 sanm.h sanm_inner.c -DDEBUG
sanm_left.o: sanm.h sanm_left.c
	cc -c -g -Wall -std=c99 sanm.h sanm_left.c -DDEBUG
sanm_right.o: sanm.h sanm_right.c
	cc -c -g -Wall -std=c99 sanm.h sanm_right.c -DDEBUG
cal_bgk.o: cal_bgk.h cal_bgk.c
	cc -c -g -Wall -std=c99 cal_bgk.h cal_bgk.c -DDEBUG
tnsol.o: tnsol.h tnsol.c
	cc -c -g -Wall -std=c99 tnsol.h tnsol.c -DDEBUG
clean:
	rm *.o
	rm *.gch
	rm *.a
