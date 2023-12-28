include ../FindPython.mk

# Compiler with debugging flags
FC = gfortran
FFLAGS = -g -fPIC -O3 -march=native -mtune=native

all: simple_cyl.$(TARGET_SUFFIX) clean_objects

clean_objects:
	rm -f *.o

simple_cyl.$(TARGET_SUFFIX): simple_cyl.f90
	LDFLAGS=-Wl,-rpath,. f2py -c -m simple_cyl simple_cyl.f90 -L. -lstell

clean:
	rm -f *.so *.o *.mod
