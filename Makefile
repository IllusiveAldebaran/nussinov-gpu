TARGET := nussinov
SRC := main.cpp nussinov.cu
HEADER := nussinov.cuh
INPUT_SEQ = inputs/ec16s.seq
FILELIST = file-list.txt
SEQ_N ?= 30


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
