FC = gfortran
FFLAGS = -g -fPIC -O3 -fopenmp -march=native -mtune=native

FFLAGS += -Wall -Wuninitialized -Wno-maybe-uninitialized -Wno-unused-label -Wno-unused-dummy-argument -Werror -Wfatal-errors -I$(HOME)/bin/libstell_dir

all: letter-canonical.x test.x test_large.x test_biotsavart.x

letter-canonical.x: canonical.o magfie.o main.f90
	$(FC) $(FFLAGS) -o $@ $^ -L. -lstell

test.x: canonical.o magfie.o test_util.f90 test.f90
	$(FC) $(FFLAGS) -o $@ $^ -L. -lstell

test_large.x: canonical.o magfie.o test_util.f90 test_large.f90
	$(FC) $(FFLAGS) -o $@ $^ -L. -lstell

test_biotsavart.x: test_biotsavart.f90 test_util.f90 biotsavart.o
	$(FC) $(FFLAGS) -o $@ $^ -L. -lstell

biotsavart.o: biotsavart.f90
	$(FC) $(FFLAGS) -c $^

canonical.o: canonical.f90 magfie.o
	$(FC) $(FFLAGS) -c $^

magfie.o: magfie_stell.f90
	$(FC) $(FFLAGS) -c $^

clean_objects:
	rm -f *.o

clean:
	rm -f *.x *.so *.o *.mod
