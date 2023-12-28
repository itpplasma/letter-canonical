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
