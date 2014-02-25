measurements.o: 	measurements.c measurements.h
	gcc -O3 -c -fpic -o $@ measurements.c

measurements_wrap.o:	measurements_wrap.c
	gcc -O3 -c -fpic -o $@ measurements_wrap.c

measurements_wrap.c:	measurements.c measurements.h measurements.i
	swig -tcl8 measurements.i

chapeau.o:	chapeau.c chapeau.h
	gcc -O3 -c -fpic -o $@ chapeau.c

cfacv.so:  cfacv.o cfacv_wrap.o measurements.o centers.o chapeau.o
	gcc -O3 -shared -o $@ cfacv.o cfacv_wrap.o measurements.o centers.o chapeau.o -lm -lgsl -lgslcblas

cfacv.o:      cfacv.c cfacv.h
	gcc -O3 -c -fpic -o $@ cfacv.c

cfacv_wrap.o: cfacv_wrap.c
	gcc -O3 -c -fpic -o $@ $?

cfacv_wrap.c: cfacv.c cfacv.i cfacv.h
	swig -tcl8 cfacv.i

centers.o:      centers.c centers.h
	gcc -O3 -c -fpic -o $@ centers.c

wrapcoords.o:	wrapcoords.c wrapcoords.h
	gcc -O3 -c -fpic -o $@ wrapcoords.c

wrapcoords_wrap.o:	wrapcoords_wrap.c
	gcc -O3 -c -fpic -o $@ $?

wrapcoords_wrap.c:	wrapcoords.c wrapcoords.h wrapcoords.i
	swig -tcl8 wrapcoords.i

clean:
	rm *.o *.so *~
