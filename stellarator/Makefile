include ../FindPython.mk

# Compiler with debugging flags
FC = gfortran
FFLAGS = -g -fPIC -O3 -march=native -mtune=native

all: simple_cyl.$(TARGET_SUFFIX) clean_objects coils.c09r00

clean_objects:
	rm -f *.o

simple_cyl.$(TARGET_SUFFIX): simple_cyl.f90
	LDFLAGS=-Wl,-rpath,. f2py -c -m simple_cyl simple_cyl.f90 -L. -lstell -lblas -llapack -I$(HOME)/bin/libstell_dir

coils.c09r00:
	wget https://princetonuniversity.github.io/STELLOPT/examples/coils.c09r00

clean:
	rm -f simple_cyl.$(TARGET_SUFFIX) *.o *.mod
