.SECONDARY:

# This makefile will compile all c stuff in a shared object
# that, thanks to SWIG, is ready to be loaded in tcl scripts

ifeq ($(CC),icc)
  LIBS+=-mkl
  # LIBS+=-lmkl_avx
  
  # LIBS+=-lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
endif

ifeq ($(CC),gcc)
   
  # CFLAGS+=--debug
  # CFLAGS+=-ggdb3
  # CFLAGS+=-O0 
  # CFLAGS+=-O3 
# CFLAGS+=-static
# CFLAGS+=-Wall
# CFLAGS+=-Wpedantic
#
# Profile
# CFLAGS+=-Wall -g -rdynamic
#
# LIBS += -lgsl -lgslcblas
LIBS += -llapack
LIBS += -lm
endif

LIBS+=-ltcl8.6

OBJ += measurements.o
OBJ += chapeau.o
OBJ += centers.o
OBJ += cfacv.o
OBJ += cvs.o

GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always 2>/dev/null || cat version 2>/dev/null || echo unknown )
CFLAGS+=-DVERSION=\"$(GIT_VERSION)\"


all: cfacv.so cvs.so chapadd chapsim catbinsp
 
%.so:  $(OBJ) %_wrap.o
	$(CC) $(CFLAGS) -shared -o $@ $^ $(LIBS)
        
%.o: %.c
	$(CC) -c -fpic $(CFLAGS) $< $(LIBS)
 
%_wrap.c: %.i %.c %.h
	swig -tcl8 $<
 
version: force
	@echo '$(GIT_VERSION)' | cmp -s - $@ || echo '$(GIT_VERSION)' > $@


catbinsp: catbinsp.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
 
chapsim: simetrize.o chapeau.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

chapadd: addchapeaus.o chapeau.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

measurements.o : measurements.h
chapeau.o : chapeau.h
simetrize.o : chapeau.h
addchapeaus.o : chapeau.h
centers.o : centers.h
cfacv.o : cfacv.h version
cvs.o : measurements.o cvs.h chapeau.h
        

.PHONY: clean force

clean:
	rm -f *.o *.so *_wrap.c
