BUILD_DIR := build
CMAKE_CACHE := $(BUILD_DIR)/CMakeCache.txt
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all
all: ninja

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(CMAKE_CACHE): CMakeLists.txt | $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake -G Ninja ..
	touch $(CMAKE_CACHE)

.PHONY: ninja
ninja: $(BUILD_NINJA)
	cd $(BUILD_DIR) && ninja

$(BUILD_NINJA): $(CMAKE_CACHE)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)
