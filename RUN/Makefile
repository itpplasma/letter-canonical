export LD_LIBRARY_PATH := $(LD_LIBRARY_PATH):..

all: Aph_of_xc.out orbit.out orbit_rmp.out

MAIN = ../BUILD/letter_canonical.x

orbit_rmp.out: letter_canonical_rmp.in $(MAIN)
	$(MAIN) letter_canonical_rmp.in

expl_impl_euler.out: field_divB0.inp letter_canonical.in ../letter_canonical.x


Aph_of_xc.out: field_divB0.inp ../test_magfie.x
	../build/test_magfie.x

clean:
	rm -f *.out fort.*
