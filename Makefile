TARGET := nussinov
SRC := main.cpp nussinov.cu
HEADER := nussinov.cuh
INPUT_SEQ = inputs/ec16s.seq
FILELIST = file-list.txt
SEQ_N ?= 30

# Default specs
BLOCK_SIZE ?= 64
MIN_LOOP_LENGTH ?= 4
GRID_SIZE ?= 1

# Add our default options here
CXXFLAGS += -DBLOCK_SIZE -DMIN_LOOP_LENGTH -DGRID_SIZE

# == compiler and flags
CXX := nvcc

CXXFLAGS := -O3 -std=c++20 -DSEQ_N=$(SEQ_N)
LDFLAGS  :=

ifeq ($(DEBUG),1)
CXXFLAGS += -DDEBUG
endif

ifeq ($(CPU),1)
CXXFLAGS += -DCPU_TARGET
endif

# == Build Rules ==

all: $(TARGET)

#$(shell head -c ${SEQ_N} ${INPUT_SEQ})

run: $(TARGET) $(HEADER)
	./$(TARGET) $(FILELIST)
	
$(TARGET): #$(OBJ)
	$(CXX) $(CXXFLAGS) $(SRC) -o $@ $(LDFLAGS)


clean:
	rm -f $(TARGET)

.PHONY: all run clean
