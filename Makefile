# Compiler with debugging flags
FC = gfortran
FFLAGS = -g -fPIC

# List of source files
SOURCES = field_divB0.f90 bdivfree.f90 tetra_grid_settings_mod.f90

# Prepend ../GORILLA/SRC/ to each source file
SOURCES := $(addprefix ../GORILLA/SRC/,$(SOURCES))

all: libfield_gorilla.so

# Target for shared library containing GORILLA magnetic field routines
libfield_gorilla.so: $(SOURCES)
	f2py -c -m field_gorilla $(SOURCES) --fcompiler=$(FC) --f90flags="$(FFLAGS)" --f77flags="$(FFLAGS)"

clean:
	rm -f *.so *.o *.mod
