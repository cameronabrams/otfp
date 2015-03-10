# This makefile will compile all c stuff in a shared object
# that, thanks to SWIG, is ready to be loaded in tcl scripts

CC=gcc

ifeq ($(CC),icc)
CCFLAGS+=-static-libgcc
LIBS += -L/usr/lib64
endif

ifeq ($(CC),gcc)
CCFLAGS+=-O3 
# CCFLAGS+=-static
# CCFLAGS=-g -Wall -Wpedantic
#
# Profile
# CCFLAGS+=-Wall -g -rdynamic
endif

# On stampede
# LIBS += -I$(TACC_GSL_INC) -L$(TACC_GSL_LIB)


LIBS += -lgsl -lgslcblas
LIBS += -lm


OBJ += measurements.o
OBJ += chapeau.o
OBJ += centers.o
OBJ += cfacv.o
OBJ += cvs.o
OBJ += cfacv_wrap.o

GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always 2>/dev/null || cat version 2>/dev/null || echo unknown )
CCFLAGS+=-DVERSION=\"$(GIT_VERSION)\"


all: cfacv.so chapadd catbinsp

cfacv.so:  $(OBJ)
	$(CC) $(CCFLAGS) -shared -o $@ $^ $(LIBS)
        
%.o: %.c
	$(CC) -c -fpic $(CCFLAGS) $< $(LIBS)
 
cfacv_wrap.c: cfacv.c cfacv.i cfacv.h
	swig -tcl8 cfacv.i

version: force
	@echo '$(GIT_VERSION)' | cmp -s - $@ || echo '$(GIT_VERSION)' > $@


catbinsp: catbinsp.o
	$(CC) $(CCFLAGS) -o $@ $^ $(LIBS)
 
chapadd: addchapeaus.o chapeau.o
	$(CC) $(CCFLAGS) -o $@ $^ $(LIBS)

measurements.o : measurements.h
chapeau.o : chapeau.h
addchapeaus.o : chapeau.h
centers.o : centers.h
cfacv.h : cvs.h measurements.h chapeau.h centers.h
cfacv.o : cfacv.h version
cvs.h : measurements.h
cvs.o : cvs.h cfacv.h 
        

.PHONY: clean force

clean:
	rm -f *.o *.so cfacv_wrap.c
