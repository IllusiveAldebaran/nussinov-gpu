
#include <string>

#define MIN_LOOP_LENGTH 4

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 256
#endif

#ifndef GRID_SIZE
#define GRID_SIZE 1
#endif


#ifdef __CUDACC__
#define HD __host__ __device__ __forceinline__
#else
#define HD inline
#endif

// Convert to triangular indexing
HD int triInd(int i, int j, int N) {
  int index = (N-MIN_LOOP_LENGTH-1)*i-((i)*(i+1))/2+j-MIN_LOOP_LENGTH-1;

  return index;
}

#ifdef __CUDACC__
#define HD __host__ __device__ __forceinline__
#else
#define HD inline
#endif
HD bool pair_check(const uint8_t* seq, int i, int j) {
  // retrieve 2 bits from a byte with mask
  uint8_t nuc1 = (seq[i / 4] >> ((i % 4) * 2)) & 0x03; 
  uint8_t nuc2 = (seq[j / 4] >> ((j % 4) * 2)) & 0x03;
  return (nuc1 ^ nuc2) == 3;
}


// Pack the string into a 2-bit array
// While 'seq' is singular it is stored contiguously in implementation
HD void encSeq(std::string unencoded_seq, uint8_t* seq, uint32_t N){

#ifdef DEBUG
  std::cout << "Length: " << N << " Seq: " << unencoded_seq << "\n" << std::endl;
#endif

  for(int i = 0; i < N; i++) {
    int8_t val = -1;
    switch(unencoded_seq[i]) { // A is 00
        case 'A': case 'a': val = 0; break; // 00
        case 'C': case 'c': val = 1; break; // 01
        case 'G': case 'g': val = 2; break; // 10
        case 'T': case 't': val = 3; break; // 11
        case 'U': case 'u': val = 3; break; // 11
    }
    seq[i / 4] |= (val << ((i % 4) * 2)); // pack 2 bits into a byte
  }
}

void nussinov_gpu_wrap(uint8_t* seq, int* DP, int N);
