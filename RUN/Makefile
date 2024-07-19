export LD_LIBRARY_PATH := $(LD_LIBRARY_PATH):..

all: Aph_of_xc.out orbit.out orbit_rmp.out

orbit_rmp.out: letter_canonical_rmp.in ../letter_canonical.x
	../letter_canonical.x letter_canonical_rmp.in

orbit.out: letter_canonical.in ../letter_canonical.x
	../letter_canonical.x

Aph_of_xc.out: field_divB0.inp ../test_magfie.x
	../test_magfie.x

clean:
	rm -f *.out fort.*
