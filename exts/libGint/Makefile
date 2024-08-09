# Makefile for C++ Program

# Compiler
CXX := nvc++

# Interpreter
PY := python3

# Compiler flags
#PROFLAGS := -p -pg
ERRFLAGS := -Wall -Wextra -Werror -Wshadow -Wformat
INCFLAGS := -I /usr/include/mkl/
LIBFLAGS := -lcudart -lcublas -lblas -lmkl_rt
CXXFLAGS := -fPIC -O3 -tp=native -std=c++14 -gpu=ccnative  -cuda $(PROFLAGS) $(ERRFLAGS) $(INCFLAGS) $(LIBFLAGS)

# Directories
SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin
SCP_DIR := tests

# Source files
CXX_SRC := $(wildcard $(SRC_DIR)/*.cpp)
CU_SRC := $(wildcard $(SRC_DIR)/*.cu)

# Object files
CXX_OBJ := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(CXX_SRC))
CU_OBJ := $(patsubst $(SRC_DIR)/%.cu,$(OBJ_DIR)/%.o,$(CU_SRC))

# Target executable
TARGET := $(BIN_DIR)/ERI_TEST

# Reference file generated using pyscf in script/CreateReference.py
REFERENCE_FILE := $(SCP_DIR)/reference.dat
PERFORMANCE_FILE := $(SCP_DIR)/performance_test.dat

# Rule to compile C++ source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c -x cu $< -o $@

# Rule to compile CUDA source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) --device-c -c $< -o $@

# Rule to link object files into the executable
$(TARGET): $(CXX_OBJ) $(CU_OBJ)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $^ -o $@

reference:
	$(PY) $(SCP_DIR)/CreateReference.py ;

# Create reference file if needed, then run test against reference values
perf: $(TARGET)
	@if [ ! -f $(PERFORMANCE_FILE) ]; then \
		$(PY) $(SCP_DIR)/PreparePerformance.py ; \
	fi
	$(TARGET) < $(PERFORMANCE_FILE)

# Create reference file if needed, then run test against reference values
test: $(TARGET)
	@if [ ! -f $(REFERENCE_FILE) ]; then \
		$(PY) $(SCP_DIR)/CreateReference.py ; \
	fi
	$(TARGET) < $(REFERENCE_FILE)

#Profile reference
profile_perf: $(TARGET)
	nsys profile --trace=cuda,nvtx $(TARGET) < $(PERFORMANCE_FILE)

profile: profile_perf

#Profile reference
profile_test: $(TARGET)
	nsys profile --trace=cuda,nvtx $(TARGET) < $(REFERENCE_FILE)

# Default rule
all: $(TARGET)

# Rule to clean the object files and the executable
clean:
	$(RM) $(CXX_OBJ) $(CU_OBJ) $(TARGET)

.PHONY: all clean test perf profile_perf profile_test

