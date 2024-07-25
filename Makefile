BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all
all: $(BUILD_NINJA)
	@$(MAKE) ninja

$(BUILD_NINJA): $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake -G Ninja .. || true

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

.PHONY: ninja
ninja: $(BUILD_NINJA)
	cd $(BUILD_DIR) && ninja

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)
