CXX := clang++
AR  := ar

CXXFLAGS := -Wall -Wextra -std=c++20 -O2 -fPIC

INC := \
	include      \
	include/core \
	include/math

SRC_DIR   := src
BUILD_DIR := build
LIB       := $(BUILD_DIR)/liblinal.a

SRC := \
    src/core/bigint/bigint.cpp \
    src/core/bigint/arithmetics.cpp \
    src/core/bigint/multiplication.cpp \
    src/core/bigint/division.cpp \
    src/core/bigint/bitwise.cpp \
    src/core/bigint/logical.cpp \
    src/core/bigint/helpers.cpp \
    src/core/bigint/io.cpp \
    src/core/bigint/advanced.cpp \
    src/core/bigint/unsafe.cpp \
    src/core/bigfloat/bigfloat.cpp \
    src/core/bigfloat/arithmetics.cpp \
    src/core/bigfloat/exponential.cpp \
    src/core/bigfloat/logarithmic.cpp \
    src/core/bigfloat/trigonometry.cpp \
    src/core/bigfloat/logical.cpp \
    src/core/bigfloat/helpers.cpp \
    src/core/bigfloat/io.cpp \
    src/math/vector.cpp \
    src/math/matrix.cpp 
    # src/math/line_2d.cpp \
    # src/math/line_nd.cpp

OBJ := $(SRC:%.cpp=$(BUILD_DIR)/%.o)

.PHONY: all clean

all: $(LIB)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(addprefix -I,$(INC)) -c $< -o $@

$(LIB): $(OBJ)
	@echo "Archiving liblinal.a"
	@$(AR) rcs $@ $^

clean:
	rm -rf $(BUILD_DIR)
