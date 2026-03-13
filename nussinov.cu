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
          if (pair_check(seq, t, j)) {
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

__global__ void traceback_gpu(cell_ind* structure, const int* DP_tri, const uint8_t* seq, int* trace_len, int N, int* stack_i, int* stack_j) {
  // sequential process- simulating recursion with stack
  if (threadIdx.x != 0 || blockIdx.x != 0) return; 

  int top = 0;

  stack_i[top] = 0;
  stack_j[top] = N - 1;
  top++;

  while (top > 0) {
    // Pop
    top--;
    int i = stack_i[top];
    int j = stack_j[top];
    if (j <= i) continue;

    if (triInd_safe(i, j, N, DP_tri) == triInd_safe(i, j - 1, N, DP_tri)) {
      stack_i[top] = i;
      stack_j[top] = j - 1;
      top++;
    } else {
      for (int k = i; k < j - MIN_LOOP_LENGTH; k++) {
        if (pair_check(seq, k, j)) {
          if (k - 1 < 0) {
            if (triInd_safe(i, j, N, DP_tri) == triInd_safe(k + 1, j - 1, N, DP_tri) + 1) {
              int idx = *trace_len;
              *trace_len = idx + 1;
              structure[idx].x = k;
              structure[idx].y = j;
                      
              // Push sub-problem
              stack_i[top] = k + 1;
              stack_j[top] = j - 1;
              top++;
              break;
            }
          } else if (triInd_safe(i, j, N, DP_tri) == (triInd_safe(i, k - 1, N, DP_tri) + triInd_safe(k + 1, j - 1, N, DP_tri) + 1)) {
            int idx = *trace_len;
            *trace_len = idx + 1;
            structure[idx].x = k;
            structure[idx].y = j;
                
            // right side
            stack_i[top] = k + 1;
            stack_j[top] = j - 1;
            top++;
                
            // left side
            stack_i[top] = i;
            stack_j[top] = k - 1;
            top++;
            break;
          }
        }
      }
    }
  } 
}


void nussinov_gpu_wrap(uint8_t* seq, cell_ind* structure, int* trace_len, int N) {
  // cudaMalloc() ?
  // DP calculation?
  uint8_t* d_seq;
  int* d_DP;

  // int* h_DP_upT = (int*)malloc(( ((N-MIN_LOOP_LENGTH)*(N-MIN_LOOP_LENGTH-1)) /2 )* sizeof(int)); // upper triangular seq

  // Allocate for DP matrix and Sequence
  CUDA_CHECK(cudaMalloc(&d_seq, N * sizeof(char)));
  CUDA_CHECK(cudaMalloc(&d_DP, ( ((N-MIN_LOOP_LENGTH)*(N-MIN_LOOP_LENGTH-1)) /2 )* sizeof(int)));

  // Copy DP matrix and Sequence
  CUDA_CHECK(cudaMemcpy(d_seq, seq, N * sizeof(char), cudaMemcpyHostToDevice));

  nussinov_gpu<<<GRID_SIZE, BLOCK_SIZE>>>(d_seq, d_DP, N);
  CUDA_CHECK(cudaDeviceSynchronize());

  // CUDA_CHECK(cudaMemcpy(h_DP_upT, d_DP, ( ((N-MIN_LOOP_LENGTH)*(N-MIN_LOOP_LENGTH-1)) /2 ) * sizeof(int), cudaMemcpyDeviceToHost));
  
  cell_ind* d_structure;
  int* d_trace_len;
  int* d_stack_i;
  int* d_stack_j;

  CUDA_CHECK(cudaMalloc(&d_structure, 2 * N * sizeof(cell_ind)));
  CUDA_CHECK(cudaMalloc(&d_trace_len, sizeof(int)));
  CUDA_CHECK(cudaMemset(d_trace_len, 0, sizeof(int)));

  CUDA_CHECK(cudaMalloc(&d_stack_i, N * sizeof(int))); 
  CUDA_CHECK(cudaMalloc(&d_stack_j, N * sizeof(int)));

  traceback_gpu<<<1, 1>>>(d_structure, d_DP, d_seq, d_trace_len, N, d_stack_i, d_stack_j);
  CUDA_CHECK(cudaDeviceSynchronize());

  CUDA_CHECK(cudaMemcpy(trace_len, d_trace_len, sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(structure, d_structure, (*trace_len) * sizeof(cell_ind), cudaMemcpyDeviceToHost));

  CUDA_CHECK(cudaFree(d_seq));
  CUDA_CHECK(cudaFree(d_DP));
  CUDA_CHECK(cudaFree(d_structure));
  CUDA_CHECK(cudaFree(d_trace_len));
  CUDA_CHECK(cudaFree(d_stack_i));
  CUDA_CHECK(cudaFree(d_stack_j));
}
