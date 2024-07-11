FC = gfortran
#FFLAGS = -g -fPIC -O0
FFLAGS = -g -fPIC -O3 -fopenmp -march=native -mtune=native -L. -Wl,-rpath=.

FFLAGS += -L. -Wl,-rpath=.
FFLAGS += -Wall -Wuninitialized -Wno-maybe-uninitialized -Wno-unused-label \
	-Wno-unused-dummy-argument -fmax-errors=1 \
	-I$(CODE)/libneo/build -L$(CODE)/libneo/build

SOURCES := magfie.f90 \
	magfie_test.f90 \
	magfie_tok.f90 \
	magfie_factory.f90 \
	canonical.f90 \
	field_can_base.f90 \
	field_can_test.f90 \
	field_can_cyl.f90 \
	field_can.f90 \
	integrator_base.f90 \
	integrator_rk45.f90 \
	integrator_euler0.f90 \
	integrator_euler1.f90 \
	integrator.f90

OBJECTS := $(SOURCES:.f90=.o)
OBJECTS := $(addprefix OBJS/, $(OBJECTS))

LIBNEO_SOURCES := libneo_kinds.f90 math_constants.f90 spl_three_to_five.f90 \
	odeint_rkf45.f90 contrib/rkf45.f90 interpolate.f90
LIBNEO_SOURCES := $(addprefix ../libneo/src/, $(LIBNEO_SOURCES))

MAGFIE_SOURCES := spline5_RZ.f90 \
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
MAGFIE_SOURCES := $(addprefix ../libneo/src/magfie/, $(MAGFIE_SOURCES))

all: letter_canonical.x test_integrator.x \
	test_field_can.x test_magfie.x test.x test_large.x test_biotsavart.x

letter_canonical.x: libfield.so libcanonical.so main.f90
	$(FC) $(FFLAGS) -o $@ $^ -lfield -lcanonical

test_integrator.x: libfield.so libcanonical.so test_integrator.f90
	$(FC) $(FFLAGS) -o $@ $^ -lfield -lcanonical

test_field_can.x: libfield.so libcanonical.so test_field_can.f90
	$(FC) $(FFLAGS) -o $@ $^ -lfield -lcanonical

test_magfie.x: libfield.so libcanonical.so test_magfie.f90
	$(FC) $(FFLAGS) -o $@ $^ -lfield -lcanonical

test.x: libfield.so libcanonical.so test_util.f90 test.f90
	$(FC) $(FFLAGS) -o $@ $^ -lfield -lcanonical

test_large.x: libfield.so libcanonical.so test_util.f90 test_large.f90
	$(FC) $(FFLAGS) -o $@ $^ -lfield -lcanonical

test_biotsavart.x: libfield.so libcanonical.so test_biotsavart.f90 test_util.f90 biotsavart.o
	$(FC) $(FFLAGS) -o $@ $^ -lfield -lcanonical

biotsavart.o: biotsavart.f90
	$(FC) $(FFLAGS) -o $@ -c $^

libcanonical.so: libfield.so $(OBJECTS)
	$(FC) $(FFLAGS) -shared -o $@ $^ -lfield

OBJS/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

libfield.so: $(LIBNEO_SOURCES) $(MAGFIE_SOURCES)
	$(FC) $(FFLAGS) -shared -o $@ $^

clean_objects:
	rm -f OBJS/*.o *.o

clean:
	rm -f *.x *.so OBJS/*.o *.o *.mod
