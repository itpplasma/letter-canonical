# Compiler with debugging flags
FC = gfortran
FFLAGS = -g -fPIC -O3 -march=native -mtune=native

# Which field variant: local, libneo, GORILLA
FIELD_VARIANT ?= libneo

# Determine the Python version dynamically
PYTHON_VERSION := $(shell python3 -c "import sys; print('{}{}'.format(sys.version_info.major, sys.version_info.minor))")

# Determine the system and architecture dynamically
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# Define a variable for the target suffix from system, architecture, and version
ifeq ($(UNAME_S),Linux)
    ifeq ($(UNAME_M),x86_64)
        TARGET_SUFFIX := cpython-$(PYTHON_VERSION)-x86_64-linux-gnu.so
    else ifeq ($(UNAME_M),aarch64)
        TARGET_SUFFIX := cpython-$(PYTHON_VERSION)-aarch64-linux-gnu.so
    else
        $(error Unsupported architecture: $(UNAME_M))
    endif
else
    $(error Unsupported operating system: $(UNAME_S))
endif

# Print suffix
$(info Target suffix: $(TARGET_SUFFIX))

# List of source files depending on variant
ifeq ($(FIELD_VARIANT), local)
    SOURCES = spline5_RZ.f90 field_divB0.f90 bdivfree.f90
else ifeq ($(FIELD_VARIANT), libneo)
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
    SOURCES := $(addprefix ../libneo/src/, $(SOURCES))
else ifeq ($(FIELD_VARIANT), GORILLA)
    SOURCES = spline5_RZ.f90 field_divB0.f90 bdivfree.f90
    SOURCES := $(addprefix ../GORILLA/SRC/, $(SOURCES))
else
    $(error Unsupported field variant: $(FIELD_VARIANT))
endif

all: libfield.so my_little_magfie.$(TARGET_SUFFIX) clean_objects

clean_objects:
	rm -f *.o

my_little_magfie.$(TARGET_SUFFIX): libfield.so my_little_magfie.f90
	LDFLAGS=-Wl,-rpath,. f2py -c -m my_little_magfie my_little_magfie.f90 -L. -lfield


# Target for shared library containing magnetic field routines
libfield.so: $(SOURCES)
	$(FC) $(FFLAGS) -shared -o $@ $^

clean:
	rm -f *.so *.o *.mod
