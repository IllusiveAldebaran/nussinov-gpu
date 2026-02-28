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

__global__ void nussinov_gpu(uint8_t* seq, int* DP, int N){
  // Initialization?
  // DP calculation?

  // -----------------------------------------------------------
  // KERNEL CONFIGURATION
  // -----------------------------------------------------------
  int tx = threadIdx.x;
  int bs = blockDim.x;
  int bx = blockIdx.x;
  //int gs = gridDim.x;

  // we now calculate 2 at once. 
  int cell_value_l;
  int cell_value_r;

  
  if (bx == 0) {

    // two thread looping
    // Notice N-k-1, does 1 less just to check in pairs
    int k;
    for(k = MIN_LOOP_LENGTH+1; k < N; k+=1){
      // Loop through diagonals (with different threads for next cell in diagonal)
      for (int i = tx; i < N-k; i += bs){
        int j = i+k;

        // unpaired value
#if MIN_LOOP_LENGTH == 0
        // we do not want to go out of bounds
        if(j==0) cell_value_l = 0;
        else cell_value_l = DP[trInd(i, j-1, N)];
#else
        cell_value_l = DP[triInd(i, j-1, N)];
#endif 
        cell_value_r = 0;

        // iterates through possible pairs for a cell
        int t;
        int pairing1 = 0;
        int pairing2l = 0;
        int pairing2r = 0;


        for(t = i; t< j-MIN_LOOP_LENGTH; t++){
          // Check paired scores (if pairing exists)
          pairing1  = 0;
          pairing2l = 0;
          pairing2r = 0;

          // check left item
          if (pair_check(seq, t, j)) {

            if(i < (t-1)-MIN_LOOP_LENGTH) 
              pairing1 = DP[triInd(i, t-1, N)];
            if(t+1 < (j-1)-MIN_LOOP_LENGTH) 
              pairing2l = DP[triInd(t+1, j-1, N)];

            cell_value_l = max(cell_value_l, pairing1 + pairing2l + 1);
          }

          // check right item
          if (pair_check(seq, t, j+1)) {
            pairing2r = 0;

            if(i < (t-1)-MIN_LOOP_LENGTH) 
              pairing1 = DP[triInd(i, t-1, N)];
            if(t+1 < (j-1)-MIN_LOOP_LENGTH) 
              pairing2r = DP[triInd(t+1, j, N)];

            cell_value_r = max(cell_value_r, pairing1 + pairing2r + 1);
          }
        }

        pairing1 = 0;
        pairing2r = 0;

        // right most item needs one more item
        // Check paired scores (if pairing exists)
        // check right item
        if (pair_check(seq, t, j+1)) {
          if(i < (t-1)-MIN_LOOP_LENGTH) 
            pairing1 = DP[triInd(i, t-1, N)];
          if(t+1 < (j-1)-MIN_LOOP_LENGTH) 
            pairing2r = DP[triInd(t+1, j, N)];
          cell_value_r = max(cell_value_r, pairing1 + pairing2r + 1);
        }

        // Second value is still diagonal, but is just 1 more value than the left one

        cell_value_r = max(cell_value_l, cell_value_r);
        DP[triInd(i,j,N)] = cell_value_l;

        // Check if right cell is within bounds... meaning we should fill it
        if (j+1<N) {
          DP[triInd(i,j+1,N)] = cell_value_r;
        }
      }
    }
  }
  
  // traceback?
}


void nussinov_gpu_wrap(uint8_t* seq, int* DP, int N) {
  // cudaMalloc() ?
  // DP calculation?
  uint8_t* d_seq;
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
