BUILD_DIR := build
CMAKE_CACHE := $(BUILD_DIR)/CMakeCache.txt
BUILD_NINJA := $(BUILD_DIR)/build.ninja
CMAKE_BUILD_TYPE ?= release

.PHONY: all
all: ninja

.PHONY: debug release

debug: CMAKE_BUILD_TYPE=Debug
debug: reconfigure

release: CMAKE_BUILD_TYPE=Release
release: reconfigure

.PHONY: reconfigure
reconfigure: $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake -G Ninja -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) ..
	touch $(CMAKE_CACHE)
	$(MAKE) ninja

$(CMAKE_CACHE):
	cd $(BUILD_DIR) && cmake -G Ninja -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) ..
	touch $(CMAKE_CACHE)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

.PHONY: ninja
ninja: $(BUILD_NINJA)
	cd $(BUILD_DIR) && ninja

$(BUILD_NINJA): $(CMAKE_CACHE) | $(BUILD_DIR)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)
