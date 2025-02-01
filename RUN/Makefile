export LD_LIBRARY_PATH := $(LD_LIBRARY_PATH):..


MAIN = ../BUILD/letter_canonical.x
CASES = fig1_dop853_cyl fig1_rk45_cyl fig1_midpoint

all: Aph_of_xc.out $(addsuffix .out, $(CASES))

%.out: %.in $(MAIN)
	$(MAIN) $<


Aph_of_xc.out: field_divB0.inp ../BUILD/test_magfie.x
	../BUILD/test_magfie.x

clean:
	rm -f *.out fort.*
