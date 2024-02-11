include FindPython.mk

# Compiler with debugging flags
FC = gfortran
FFLAGS = -g -fPIC -O3 -fopenmp -march=native -mtune=native
# FFLAGS = -g -fPIC -Og -fopenmp

FFLAGS += -Wall -Wuninitialized -Wno-maybe-uninitialized -Wno-unused-label -Wno-unused-dummy-argument -Werror -Wfatal-errors

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

SOURCES += odeint_rkf45.f90 contrib/rkf45.f90 interpolate.f90
SOURCES := $(addprefix ../libneo/src/, $(SOURCES))

all: libfield.so my_little_magfie.$(TARGET_SUFFIX) letter-canonical.x \
	clean_objects

letter-canonical.x: canonical.o my_little_magfie.o main.f90
	$(FC) $(FFLAGS) -o $@ $^ -L. -lfield-lbspline-fortran

canonical.o: canonical.f90 my_little_magfie.o
	$(FC) $(FFLAGS) -c $^

my_little_magfie.o: my_little_magfie.f90 libfield.so
	$(FC) $(FFLAGS) -c $^

clean_objects:
	rm -f *.o

my_little_magfie.$(TARGET_SUFFIX): libfield.so my_little_magfie.f90
	LDFLAGS=-Wl,-rpath,. f2py -c -m my_little_magfie my_little_magfie.f90 -L. -lfield

libfield.so: $(SOURCES)
	$(FC) $(FFLAGS) -shared -o $@ $^

clean:
	rm -f *.x *.so *.o *.mod
