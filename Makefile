BUILD_DIR := BUILD
CMAKE_CACHE := $(BUILD_DIR)/CMakeCache.txt
BUILD_NINJA := $(BUILD_DIR)/build.ninja
CMAKE_BUILD_TYPE ?= release

.PHONY: all debug release clean paper ninja reconfigure
all: reconfigure

debug: CMAKE_BUILD_TYPE=Debug
debug: reconfigure

release: CMAKE_BUILD_TYPE=Release
release: reconfigure

reconfigure: $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake -G Ninja -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) ..
	touch $(CMAKE_CACHE)
	$(MAKE) ninja

$(CMAKE_CACHE):
	cd $(BUILD_DIR) && cmake -G Ninja -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) ..
	touch $(CMAKE_CACHE)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

ninja: $(BUILD_NINJA)
	cd $(BUILD_DIR) && ninja

$(BUILD_NINJA): $(CMAKE_CACHE) | $(BUILD_DIR)

paper:
	@echo "\
--------------------------------------------------------------------------- \n\
| Log in with username 'git' and your authentication token as a password. | \n\
| Generate your token on https://www.overleaf.com/user/settings           | \n\
---------------------------------------------------------------------------"
	@if [ ! -d "paper" ]; then \
		git clone https://git@git.overleaf.com/675740e8b9283f4c4678d63b paper; \
	else \
		cd paper; git pull; cd ..; \
	fi

clean:
	rm -rf $(BUILD_DIR)
