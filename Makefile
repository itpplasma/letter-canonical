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

letter-canonical.x: libfield.so libcanonical.so main.f90
	$(FC) $(FFLAGS) -o $@ $^ -L. -lstell -lfield -lcanonical

test.x: libfield.so libcanonical.so test_util.f90 test.f90
	$(FC) $(FFLAGS) -o $@ $^ -L. -lstell -lfield -lcanonical

test_large.x: libfield.so libcanonical.so test_util.f90 test_large.f90
	$(FC) $(FFLAGS) -o $@ $^ -L. -lstell -lfield -lcanonical

test_biotsavart.x: libfield.so libcanonical.so test_biotsavart.f90 test_util.f90 biotsavart.o
	$(FC) $(FFLAGS) -o $@ $^ -L. -lstell -lfield -lcanonical

biotsavart.o: biotsavart.f90
	$(FC) $(FFLAGS) -o $@ -c $^

libcanonical.so: magfie.f90 magfie_test.f90 magfie_factory.f90 canonical.f90
	$(FC) $(FFLAGS) -shared -o $@ $^

libfield.so: $(SOURCES)
	$(FC) $(FFLAGS) -shared -o $@ $^

clean_objects:
	rm -f *.o

clean:
	rm -f *.x *.so *.o *.mod
