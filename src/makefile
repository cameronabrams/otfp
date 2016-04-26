
# This makefile will compile all c stuff in a shared object
# that, thanks to SWIG, is ready to be loaded in tcl scripts

CC=gcc

ifeq ($(CC),icc)
CCFLAGS+=-static-libgcc
LIBS += -L/usr/lib64
endif

ifeq ($(CC),gcc)
# CCFLAGS+=-g
# CCFLAGS+=-O0 
# CCFLAGS+=-O3 
# CCFLAGS+=-static
# CCFLAGS+=-Wall
# CCFLAGS+=-Wpedantic
#
# Profile
# CCFLAGS+=-Wall -g -rdynamic
endif

# On stampede
# LIBS += -I$(TACC_GSL_INC) -L$(TACC_GSL_LIB)


# LIBS += -lgsl -lgslcblas
LIBS += -llapack
LIBS += -lm


OBJ += measurements.o
OBJ += chapeau.o
OBJ += centers.o
OBJ += cfacv.o
OBJ += cvs.o

GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always 2>/dev/null || cat version 2>/dev/null || echo unknown )
CCFLAGS+=-DVERSION=\"$(GIT_VERSION)\"


all: cfacv.so chapadd catbinsp
 
%.so:  $(OBJ) %_wrap.o
	$(CC) $(CCFLAGS) -shared -o $@ $^ $(LIBS)
        
%.o: %.c
	$(CC) -c -fpic $(CCFLAGS) $< $(LIBS)
 
%_wrap.c: %.i %.c %.h
	swig -tcl8 $<
 
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
	rm -f *.o *.so *_wrap.c
