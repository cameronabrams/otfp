# This makefile will compile all c stuff in a shared object
# that, thanks to SWIG, is ready to be loaded in tcl scripts

CC=gcc
CCFLAGS=-O3
LIBS = -lm -lgsl -lgslcblas

OBJ += measurements.o
OBJ += chapeau.o
OBJ += centers.o
OBJ += cfacv.o
OBJ += cfacv_wrap.o

cfacv.so:  $(OBJ) cfacv_wrap.o
	$(CC) $(CCFLAGS) -shared -o $@ $^ $(LIBS)
        
%.o: %.c
	$(CC) -c -fpic $(CCFLAGS) $<
 
cfacv_wrap.c: cfacv.c cfacv.i cfacv.h
	swig -tcl8 cfacv.i

measurements.o : measurements.h
chapeau.o : chapeau.h
centers.o : centers.h
cfacf.h : measurements.h chapeau.h centers.h
cfacv.o : cfacv.h


.PHONY: clean

clean:
	rm -f *.o *.so cfacv_wrap.c
