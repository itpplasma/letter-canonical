FC = gfortran
FFLAGS = -g -fPIC -O3 -fopenmp -march=native -mtune=native

FFLAGS += -Wall -Wuninitialized -Wno-maybe-uninitialized -Wno-unused-label \
	-Wno-unused-dummy-argument -Werror -Wfatal-errors -fmax-errors=1 \
	-I$(HOME)/bin/libstell_dir -L$(HOME)/bin/libstell_dir \
	-I$(CODE)/libneo/build -L$(CODE)/libneo/build

SOURCES = libneo_kinds.f90 math_constants.f90 spl_three_to_five.f90

MAGFIE_SOURCES = spline5_RZ.f90

SOURCES += $(addprefix magfie/, $(MAGFIE_SOURCES))

SOURCES += odeint_rkf45.f90 contrib/rkf45.f90 interpolate.f90
SOURCES := $(addprefix ../libneo/src/, $(SOURCES))

all: letter-canonical.x test.x test_large.x test_biotsavart.x

letter-canonical.x: libfield.so canonical.o magfie.o main.f90
	$(FC) $(FFLAGS) -o $@ $^ -L. -lstell -lfield

test.x: libfield.so canonical.o magfie.o test_util.f90 test.f90
	$(FC) $(FFLAGS) -o $@ $^ -L. -lstell -lfield

test_large.x: libfield.so canonical.o magfie.o test_util.f90 test_large.f90
	$(FC) $(FFLAGS) -o $@ $^ -L. -lstell -lfield

test_biotsavart.x: libfield.so test_biotsavart.f90 test_util.f90 biotsavart.o
	$(FC) $(FFLAGS) -o $@ $^ -L. -lstell -lfield

biotsavart.o: biotsavart.f90
	$(FC) $(FFLAGS) -o $@ -c $^

canonical.o: canonical.f90 magfie.o
	$(FC) $(FFLAGS) -o $@ -c $^

magfie.o: magfie_stell.f90
	$(FC) $(FFLAGS) -o $@ -c $^

libfield.so: $(SOURCES)
	$(FC) $(FFLAGS) -shared -o $@ $^

clean_objects:
	rm -f *.o

clean:
	rm -f *.x *.so *.o *.mod
