include ../FindPython.mk

# Compiler with debugging flags
FC = gfortran
FFLAGS = -g -fPIC -O3 -march=native -mtune=native -I$(HOME)/bin/libstell_dir

all: simple_cyl.$(TARGET_SUFFIX) clean_objects coils.c09r00 \
	letter-canonical-stel.x

clean_objects:
	rm -f *.o

letter-canonical-stel.x: canonical.o my_little_magfie.o main.f90
	$(FC) $(FFLAGS) -o $@ $^ -L. -lstell -lfield

canonical.o: ../canonical.f90 my_little_magfie.o
	$(FC) $(FFLAGS) -c $^

my_little_magfie.o: my_little_magfie.f90 libstell.so
	$(FC) $(FFLAGS) -c $^

simple_cyl.$(TARGET_SUFFIX): my_little_magfie.f90
	LDFLAGS=-Wl,-rpath,. f2py -c -m simple_cyl my_little_magfie.f90 -L$(HOME)/bin -lstell -lblas -llapack -I$(HOME)/bin/libstell_dir

coils.c09r00:
	wget https://princetonuniversity.github.io/STELLOPT/examples/coils.c09r00

clean:
	rm -f simple_cyl.$(TARGET_SUFFIX) *.o *.mod
