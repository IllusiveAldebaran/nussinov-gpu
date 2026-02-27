TARGET := nussinov
SRC := main.cpp nussinov.cu
INPUT_SEQ = inputs/ec16s.seq
SEQ_N ?= 30


# == compiler and flags

CXX := nvcc

CXXFLAGS := -O3 -std=c++20 
LDFLAGS  :=

ifeq ($(DEBUG),1)
CXXFLAGS += -DDEBUG
endif

ifeq ($(CPU),1)
CXXFLAGS += -DCPU_TARGET
endif

# == Build Rules ==

all: run

run: $(TARGET)
	./$(TARGET) $(shell head -c ${SEQ_N} ${INPUT_SEQ})
	
$(TARGET): #$(OBJ)
	$(CXX) $(CXXFLAGS) $(SRC) -o $@ $(LDFLAGS)


clean:
	rm -f $(TARGET)

.PHONY: all run clean
