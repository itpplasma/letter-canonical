FC = gfortran
FFLAGS = -g -fPIC -O3 -fopenmp -march=native -mtune=native -L.

FFLAGS += -Wall -Wuninitialized -Wno-maybe-uninitialized -Wno-unused-label \
	-Wno-unused-dummy-argument -Werror -Wfatal-errors -fmax-errors=1 \
	-I$(HOME)/bin/libstell_dir -L$(HOME)/bin/libstell_dir \
	-I$(CODE)/libneo/build -L$(CODE)/libneo/build \
	-I$(CODE)/external/bspline-fortran/build/lib \
	-L$(CODE)/external/bspline-fortran/build/lib

SOURCES = libneo_kinds.f90 math_constants.f90 spl_three_to_five.f90

MAGFIE_SOURCES = spline5_RZ.f90 \
theta_rz_mod.f90 \
extract_fluxcoord_mod.f90 \
amn_mod.f90 \
field_mod.f90 \
field_eq_mod.f90 \
field_c_mod.f90 \
input_files.f90 \
input_files.f90 \
inthecore_mod.f90 \
field_divB0.f90 \
bdivfree_mod.f90 \
bdivfree.f90

SOURCES += $(addprefix magfie/, $(MAGFIE_SOURCES))

SOURCES += odeint_rkf45.f90 contrib/rkf45.f90
SOURCES := $(addprefix ../libneo/src/, $(SOURCES))

all: letter-canonical.x test.x test_large.x test_biotsavart.x

letter-canonical.x: libfield.so libcanonical.so main.f90
	$(FC) $(FFLAGS) -o $@ $^ -lstell -lfield -lcanonical -lbspline-fortran

test.x: libfield.so libcanonical.so test_util.f90 test.f90
	$(FC) $(FFLAGS) -o $@ $^ -lstell -lfield -lcanonical -lbspline-fortran

test_large.x: libfield.so libcanonical.so test_util.f90 test_large.f90
	$(FC) $(FFLAGS) -o $@ $^ -lstell -lfield -lcanonical -lbspline-fortran

test_biotsavart.x: libfield.so libcanonical.so test_biotsavart.f90 test_util.f90 biotsavart.o
	$(FC) $(FFLAGS) -o $@ $^ -lstell -lfield -lcanonical -lbspline-fortran

biotsavart.o: biotsavart.f90
	$(FC) $(FFLAGS) -o $@ -c $^

libcanonical.so: interpolate_adapter.f90 magfie.f90 magfie_test.f90 magfie_tok.f90 magfie_factory.f90 canonical.f90 libfield.so
	$(FC) $(FFLAGS) -shared -o $@ $^ -lfield -lbspline-fortran

libfield.so: $(SOURCES)
	$(FC) $(FFLAGS) -shared -o $@ $^

clean_objects:
	rm -f *.o

clean:
	rm -f *.x *.so *.o *.mod
