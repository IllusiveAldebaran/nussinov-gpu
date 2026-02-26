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
__global__ void nussinov_gpu(char* seq, int* DP, int N){
  // Initialization?
  // DP calculation?

  // -----------------------------------------------------------
  // KERNEL CONFIGURATION
  // -----------------------------------------------------------
  int tx = threadIdx.x;
  int bs = blockDim.x;
  int bx = blockIdx.x;
  int gs = gridDim.x;

  int cell_value;


  
  if (bx == 0) {
    // initialization
    // NxN matrix with scores of optimal pairings
    for(int k =0; k< MIN_LOOP_LENGTH; k++) {
      int j;
      for(int i = 0; i < N-k; i++){
        j = i + k;
        DP[N*i+j] = 0; //INT_MIN;
      }

    }



    for(int k = MIN_LOOP_LENGTH; k < N; k++){
      // TODO: Per Thread Diagonal Looping here
      for (int i = tx; i < N-k; i += bs){
        int j = i+k;


        cell_value = 0;
        if (i < j-MIN_LOOP_LENGTH) {

          // unpaired value
#if MIN_LOOP_LENGTH == 0
          // we do not want to go out of bounds
          if(j==0) cell_value = 0;
          else cell_value = DP[N*i + (j-1)];
#else
          cell_value = DP[N*i + (j-1)];
#endif 

          // iterates through possible pairs for a cell
          for(int t = i; t< j-MIN_LOOP_LENGTH; t++){
            // Check paired scores (if pairing exists)
            if (pair_check_gpu(seq[t], seq[j])) {
              int pairing1 = 0;
              int pairing2 = 0;

              if(i < (t-1)-MIN_LOOP_LENGTH) 
                pairing1 = DP[N*i + (t-1)];
              if(t+1 < (j-1)-MIN_LOOP_LENGTH) 
                pairing2 = DP[N*(t+1) + (j-1)];

              cell_value = max(cell_value, pairing1 + pairing2 + 1);
            }
          }
        }

        DP[N*i+j] = cell_value;
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

  // Allocate for DP matrix and Sequence
  CUDA_CHECK(cudaMalloc(&d_seq, N * sizeof(char)));
  CUDA_CHECK(cudaMalloc(&d_DP, N * N * sizeof(int)));

  // Copy DP matrix and Sequence
  CUDA_CHECK(cudaMemcpy(d_seq, seq, N * sizeof(char), cudaMemcpyHostToDevice));

  nussinov_gpu<<<GRID_SIZE, BLOCK_SIZE>>>(d_seq, d_DP, N);

  CUDA_CHECK(cudaMemcpy(DP, d_DP, N * N * sizeof(int), cudaMemcpyDeviceToHost));

  CUDA_CHECK(cudaDeviceSynchronize());

  CUDA_CHECK(cudaFree(d_seq));
  CUDA_CHECK(cudaFree(d_DP));

  // traceback?
  // cudaFree() ?
}
