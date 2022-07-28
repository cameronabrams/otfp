.SECONDARY:

# This makefile will compile all c stuff in a shared object
# that, thanks to SWIG, is ready to be loaded in tcl scripts

ifeq ($(CC),icc)
  LIBS+=-mkl
  # LIBS+=-lmkl_avx
  
  # LIBS+=-lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
endif
   
  CFLAGS+=--debug
  CFLAGS+=-g
  CFLAGS+=-O0 
  # CFLAGS+=-O3 
# CFLAGS+=-static
# CFLAGS+=-Wall
# CFLAGS+=-Wpedantic
#  
ifeq ($(CC),gcc)
  
# Profile
# CFLAGS+=-Wall -g -rdynamic
#
# LIBS += -lgsl -lgslcblas
LIBS += -llapack -lblas
LIBS += -lm
endif

LIBS+=-ltcl8.6 
# Set  /home/alexis/apps/local/librsb/lib that shuld be in the LD_LIBRARY_PATH
# LIBS+= -I/home/alexis/apps/local/librsb/include/ -lrsb 
# Set  /home/alexis/apps/git/SuiteSparse/
LIBS+= -I/home/alexis/apps/git/SuiteSparse/include/ -lumfpack -lamd -lcholmod -lsuitesparseconfig -lm -llapack -lblas

OBJ += measurements.o
OBJ += chapeau.o
OBJ += centers.o
OBJ += cfacv.o
OBJ += cvs.o

GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always 2>/dev/null || cat version 2>/dev/null || echo unknown )
CFLAGS+=-DVERSION=\"$(GIT_VERSION)\"


all: cfacv.so cvs.so chapadd chapsim catbinsp chapcrop chapconv
 
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

chapconv: convchapeau.o chapeau.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

chapadd: addchapeaus.o chapeau.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

chapcrop: cropchapeaus.o chapeau.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

measurements.o : measurements.h
chapeau.o : chapeau.h chapeau_obj.h
simetrize.o : chapeau.h chapeau_obj.h
addchapeaus.o : chapeau.h chapeau_obj.h
cropchapeaus.o : chapeau.h chapeau_obj.h
convchapeaus.o : chapeau.h chapeau_obj.h chapeau_old.h
centers.o : centers.h
cfacv.o : cfacv.h version cfacv_obj.h
cvs.o : measurements.o cvs.h chapeau.h cvs_obj.h chapeau_obj.h
        

.PHONY: clean force

clean:
	rm -f *.o *.so *_wrap.c
