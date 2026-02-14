# =======================
# Toolchain
# =======================
CXX := clang++
AR  := ar

CXXFLAGS := -Wall -Wextra -std=c++20 -O2 -g
LDFLAGS  :=

SRC_DIR   := src
BUILD_DIR := build

# =======================
# Third-party library
# =======================
TP_DIR   := third_party/linal-sdk
TP_LIB   := $(TP_DIR)/build/liblinal.a

INC := include \
$(TP_DIR)/include/core \
$(TP_DIR)/include/math

# =======================
# Runtime
# =======================
DEBUGGER_CMD := pwndbg
ARGS :=

TASK ?=

ifeq ($(V),1)
	Q :=
else
	Q := @
endif

.PHONY: all clean run pwn valgrind help

# =======================
# Default
# =======================
all:
	$(Q)echo "Use: make run TASK=<n>"

# =======================
# Ensure third-party library is built
# =======================
$(TP_LIB):
	$(Q)$(MAKE) -C $(TP_DIR)

# =======================
# Build tasks
# =======================
$(BUILD_DIR)/%: $(SRC_DIR)/%.cpp $(TP_LIB)
	$(Q)echo "Building task $*"
	$(Q)mkdir -p $(BUILD_DIR)
	$(Q)$(CXX) $(CXXFLAGS) $(addprefix -I,$(INC)) $< $(TP_LIB) $(LDFLAGS) -o $@

# =======================
# Run
# =======================
run:
	@if [ -z "$(TASK)" ]; then \
		echo "Error: TASK is not set"; exit 1; \
	fi
	$(MAKE) $(BUILD_DIR)/$(TASK)
	$(Q)./$(BUILD_DIR)/$(TASK) $(ARGS)

# =======================
# Debug
# =======================
pwn:
	@if [ -z "$(TASK)" ]; then \
		echo "Error: TASK is not set"; exit 1; \
	fi
	$(MAKE) $(BUILD_DIR)/$(TASK)
	$(Q)$(DEBUGGER_CMD) ./$(BUILD_DIR)/$(TASK)

# =======================
# Valgrind
# =======================
valgrind:
	@if [ -z "$(TASK)" ]; then \
		echo "Error: TASK is not set"; exit 1; \
	fi
	$(MAKE) $(BUILD_DIR)/$(TASK)
	$(Q)valgrind --leak-check=full --show-leak-kinds=all \
		--track-origins=yes --error-exitcode=1 \
		./$(BUILD_DIR)/$(TASK) $(ARGS)

# =======================
# Clean
# =======================
clean:
	$(Q)rm -rf $(BUILD_DIR)
	$(Q)$(MAKE) -C $(TP_DIR) clean

# =======================
# Help
# =======================
help:
	$(Q)echo "Targets:"
	$(Q)echo "  run TASK=<n>       build & run src/<n>.cpp"
	$(Q)echo "  pwn TASK=<n>       debug task"
	$(Q)echo "  valgrind TASK=<n>  run under valgrind"
	$(Q)echo "  clean              remove build/"
	$(Q)echo ""
	$(Q)echo "Variables:"
	$(Q)echo "  TASK=<n>  task number"
	$(Q)echo "  ARGS=...  runtime args"
	$(Q)echo "  V=1       verbose build"
