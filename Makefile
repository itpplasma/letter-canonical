BUILD_DIR := BUILD
CMAKE_CACHE := $(BUILD_DIR)/CMakeCache.txt
BUILD_NINJA := $(BUILD_DIR)/build.ninja
CMAKE_BUILD_TYPE ?= release

.PHONY: all debug release
all: reconfigure

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

paper:
	git clone https://git@git.overleaf.com/675740e8b9283f4c4678d63b paper

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)
