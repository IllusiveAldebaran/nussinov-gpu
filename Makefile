TARGET := nussinov
SRC := main.cpp nussinov.cu
HEADER := nussinov.cuh
FILELIST = file-list.txt

# Base Compiler Flags 
CXX := nvcc
CXXFLAGS := -O3 -std=c++20

# Default specs and overrides
SEQ_N ?= 30
BLOCK_SIZE ?= 64
MIN_LOOP_LENGTH ?= 4
GRID_SIZE ?= 1

ifdef BLOCK
    BLOCK_SIZE := $(BLOCK)
endif
ifdef LOOP
    MIN_LOOP_LENGTH := $(LOOP)
endif
ifdef GRID
    GRID_SIZE := $(GRID)
endif

# Append all defines to CXXFLAGS
CXXFLAGS += -DSEQ_N=$(SEQ_N)
CXXFLAGS += -DBLOCK_SIZE=$(BLOCK_SIZE)
CXXFLAGS += -DMIN_LOOP_LENGTH=$(MIN_LOOP_LENGTH)
CXXFLAGS += -DGRID_SIZE=$(GRID_SIZE)

ifeq ($(DEBUG),1)
    CXXFLAGS += -DDEBUG
endif

ifeq ($(CPU),1)
    CXXFLAGS += -DCPU_TARGET
endif

# == Build Rules ==
all: $(TARGET)

$(TARGET): $(SRC) $(HEADER)
	$(CXX) $(CXXFLAGS) $(SRC) -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: all run clean