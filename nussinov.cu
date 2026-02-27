/*
 * List of kernels for nussinov
 * 
 */

// for the macro and structs
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "nussinov.cuh"

#define CUDA_CHECK(expression)\
{\
    const cudaError_t err = expression;\
    if (err != cudaSuccess) {\
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));\
        exit(1);\
    }\
}\

__device__ __forceinline__
bool pair_check_gpu(char nuc1, char nuc2) {
  bool check = false;

  // TODO: encode and compare via bool
  if(nuc1 == 'A' && nuc2 == 'U') check = true;
  if(nuc1 == 'U' && nuc2 == 'A') check = true;
  if(nuc1 == 'C' && nuc2 == 'G') check = true;
  if(nuc1 == 'G' && nuc2 == 'C') check = true;

  return check;
}

// Convert to triangular indexing
__host__ __device__ __forceinline__
int triInd(int i, int j, int N) {
  //int index = (N-MIN_LOOP_LENGTH-1)*i+((i)*(i+1))/2-(MIN_LOOP_LENGTH+1)*(i+1)+(j);
  int index = (N-MIN_LOOP_LENGTH-1)*i-((i)*(i+1))/2+j-MIN_LOOP_LENGTH-1;

  return index;
}


__global__ void nussinov_gpu(char* seq, int* DP, int N){
  // Initialization?
  // DP calculation?

  // -----------------------------------------------------------
  // KERNEL CONFIGURATION
  // -----------------------------------------------------------
  int tx = threadIdx.x;
  int bs = blockDim.x;
  int bx = blockIdx.x;
  //int gs = gridDim.x;

  int cell_value;

  
  if (bx == 0) {

    for(int k = MIN_LOOP_LENGTH+1; k < N; k++){
      // Loop through diagonals (with different threads for next cell in diagonal)
      for (int i = tx; i < N-k; i += bs){
        int j = i+k;

        // unpaired value
#if MIN_LOOP_LENGTH == 0
        // we do not want to go out of bounds
        if(j==0) cell_value = 0;
        else cell_value = DP[trInd(i, j-1, N)];
#else
        cell_value = DP[triInd(i, j-1, N)];
#endif 

        // iterates through possible pairs for a cell
        for(int t = i; t< j-MIN_LOOP_LENGTH; t++){
          // Check paired scores (if pairing exists)
          if (pair_check_gpu(seq[t], seq[j])) {
            int pairing1 = 0;
            int pairing2 = 0;

            if(i < (t-1)-MIN_LOOP_LENGTH) 
              pairing1 = DP[triInd(i, t-1, N)];
            if(t+1 < (j-1)-MIN_LOOP_LENGTH) 
              pairing2 = DP[triInd(t+1, j-1, N)];

            cell_value = max(cell_value, pairing1 + pairing2 + 1);
          }
        }

        DP[triInd(i,j,N)] = cell_value;
      }
    }
  }
  
  // traceback?
}


void nussinov_gpu_wrap(char* seq, int* DP, int N) {
  // cudaMalloc() ?
  // DP calculation?
  char* d_seq;
  int* d_DP;

  int* h_DP_upT = (int*)malloc(( ((N-MIN_LOOP_LENGTH)*(N-MIN_LOOP_LENGTH-1)) /2 )* sizeof(int)); // upper triangular seq

  // Allocate for DP matrix and Sequence
  CUDA_CHECK(cudaMalloc(&d_seq, N * sizeof(char)));
  CUDA_CHECK(cudaMalloc(&d_DP, ( ((N-MIN_LOOP_LENGTH)*(N-MIN_LOOP_LENGTH-1)) /2 )* sizeof(int)));

  // Copy DP matrix and Sequence
  CUDA_CHECK(cudaMemcpy(d_seq, seq, N * sizeof(char), cudaMemcpyHostToDevice));

  nussinov_gpu<<<GRID_SIZE, BLOCK_SIZE>>>(d_seq, d_DP, N);

  CUDA_CHECK(cudaMemcpy(h_DP_upT, d_DP, ( ((N-MIN_LOOP_LENGTH)*(N-MIN_LOOP_LENGTH-1)) /2 ) * sizeof(int), cudaMemcpyDeviceToHost));

  // traceback?

  // Copy uptriangular matrix to real NxN Matrix
  for(int k = MIN_LOOP_LENGTH+1; k < N; k++){
    for (int i = 0; i < N-k; i++){
      int j = i+k;
      DP[N*i+j] = h_DP_upT[triInd(i, j, N)];
    }
  }

  CUDA_CHECK(cudaDeviceSynchronize());

  CUDA_CHECK(cudaFree(d_seq));
  CUDA_CHECK(cudaFree(d_DP));
}
