/*
 * List of kernels for nussinov
 * 
 */

// for the macro and structs
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <bits/stdc++.h> // max function

#include "nussinov.cuh"

#define CUDA_CHECK(expression)\
{\
    const cudaError_t err = expression;\
    if (err != cudaSuccess) {\
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));\
        exit(1);\
    }\
}\

__global__ void nussinov_gpu(uint8_t* seqs, uint32_t* seq_offsets, uint32_t* seq_lengths, 
                             uint32_t* dp_offsets, int* batched_DP, uint32_t N){
  // -----------------------------------------------------------
  // KERNEL CONFIGURATION
  // -----------------------------------------------------------
  int tx = threadIdx.x;
  int bs = blockDim.x;
  int bx = blockIdx.x;
  //int gs = gridDim.x;

  uint32_t len = seq_lengths[bx];
  uint8_t* seq = &seqs[seq_offsets[bx]];
  int* DP = &batched_DP[dp_offsets[bx]];

  for(int k = MIN_LOOP_LENGTH+1; k < len; k++){
    // Loop through diagonals (with different threads for next cell in diagonal)
    for (int i = tx; i < len-k; i += bs){
      int j = i+k;
      int cell_value;

      if (i >= j - MIN_LOOP_LENGTH) { // prevent out of bounds
        cell_value = 0;
      } else {
        // unpaired value
#if MIN_LOOP_LENGTH == 0
        // we do not want to go out of bounds
        if(j==0) cell_value = 0;
        else cell_value = DP[triInd(i, j-1, len)];
#else
        cell_value = DP[triInd(i, j-1, len)];
#endif 

        // iterates through possible pairs for a cell
        for(int t = i; t< j-MIN_LOOP_LENGTH; t++){
          // Check paired scores (if pairing exists)
          if (pair_check(seq, t, j)) {
            int pairing1 = 0;
            int pairing2 = 0;

            if(i < (t-1)-MIN_LOOP_LENGTH) 
              pairing1 = DP[triInd(i, t-1, len)];
            if(t+1 < (j-1)-MIN_LOOP_LENGTH) 
              pairing2 = DP[triInd(t+1, j-1, len)];

            cell_value = max(cell_value, pairing1 + pairing2 + 1);
          }
        }
      }

      DP[triInd(i,j,len)] = cell_value;
    }
    __syncthreads();
  }
  
  // traceback?
}


void nussinov_gpu_wrap(uint8_t* seqs, uint32_t* seq_offsets, uint32_t* seq_lengths, uint32_t* dp_offsets, 
                        uint32_t N, uint32_t total_bytes, uint32_t total_dp_cells, int* batched_DP) {

  uint8_t* d_seqs;
  uint32_t *d_seq_offsets, *d_seq_lengths, *d_dp_offsets;
  int *d_batched_DP;

  CUDA_CHECK(cudaMalloc(&d_seqs, total_bytes * sizeof(uint8_t)));
  CUDA_CHECK(cudaMalloc(&d_seq_offsets, N * sizeof(uint32_t)));
  CUDA_CHECK(cudaMalloc(&d_seq_lengths, N * sizeof(uint32_t)));
  CUDA_CHECK(cudaMalloc(&d_dp_offsets, N * sizeof(uint32_t)));
  CUDA_CHECK(cudaMalloc(&d_batched_DP, total_dp_cells * sizeof(int)));

  CUDA_CHECK(cudaMemcpy(d_seqs, seqs, total_bytes * sizeof(uint8_t), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_seq_offsets, seq_offsets, N * sizeof(uint32_t), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_seq_lengths, seq_lengths, N * sizeof(uint32_t), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_dp_offsets, dp_offsets, N * sizeof(uint32_t), cudaMemcpyHostToDevice));
  
  CUDA_CHECK(cudaMemset(d_batched_DP, 0, total_dp_cells * sizeof(int))); //calloc

  int num_blocks = N;

  nussinov_gpu<<<num_blocks, BLOCK_SIZE>>>(d_seqs, d_seq_offsets, d_seq_lengths, d_dp_offsets, d_batched_DP, N);
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());
  CUDA_CHECK(cudaMemcpy(batched_DP, d_batched_DP, total_dp_cells * sizeof(int), cudaMemcpyDeviceToHost));

  CUDA_CHECK(cudaFree(d_seqs));
  CUDA_CHECK(cudaFree(d_seq_offsets));
  CUDA_CHECK(cudaFree(d_seq_lengths));
  CUDA_CHECK(cudaFree(d_dp_offsets));
  CUDA_CHECK(cudaFree(d_batched_DP));
}
