

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

// cell index
struct cell_ind{
  int x;
  int y;
};

// Convert to triangular indexing
HD int triInd(int i, int j, int N) {
  int index = (N-MIN_LOOP_LENGTH-1)*i-((i)*(i+1))/2+j-MIN_LOOP_LENGTH-1;

  return index;
}

HD int triInd_safe(int i, int j, int N, const int* DP) {
    if (j - i <= MIN_LOOP_LENGTH) {
        return 0; 
    }
    return DP[triInd(i, j, N)];
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

void nussinov_gpu_wrap(uint8_t* seq, cell_ind* structure, int* trace_len, int N);
