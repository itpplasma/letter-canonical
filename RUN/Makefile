export LD_LIBRARY_PATH := $(LD_LIBRARY_PATH):..


MAIN = ../BUILD/letter_canonical.x
CASES = fig1_dop853_cyl fig1_rk45_cyl fig1_midpoint \
        fig2_dop853_cyl fig2_rk45_cyl fig2_rk45_can fig2_midpoint fig2_expl_impl_euler \
		fig3_dop853_cyl fig3_rk45_cyl fig3_rk45_can fig3_midpoint fig3_expl_impl_euler

all: Aph_of_xc.out $(addsuffix .out, $(CASES))

%.out: %.in $(MAIN)
	$(MAIN) $<


Aph_of_xc.out: field_divB0.inp ../BUILD/test_magfie.x
	../BUILD/test_magfie.x

clean:
	rm -f *.out fort.*
