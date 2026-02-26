TARGET := nussinov
SRC := main.cpp # nussinov.cu
INPUT_SEQ = inputs/ec16s.seq
SEQ_N ?= 30


# == compiler and flags

CXX ?= nvcc

CXXFLAGS := -O3 -std=c++20 -Wall -Wextra
CXXFLAGS += -DFORCE_ALL_OMP -mavx2 -faligned-new -DSSE_AVX2 -DOMPGPU -fopenmp-offload-mandatory --offload-arch=native -fopenmp-force-usm
LDFLAGS  :=
#DEFINES  := -DDEBUG

# == Build Rules ==

all: run

run: $(TARGET)
	./$(TARGET) $(shell head -c ${SEQ_N} ${INPUT_SEQ})
	
$(TARGET): #$(OBJ)
	$(CXX) $(SRC) -o $@ $(LDFLAGS)

#%.o: %.cpp
#	$(CXX) $(CXXFLAGS) $(DEFINES) -c $< -o $@

clean:
	rm -f $(TARGET)

.PHONY: all run clean
