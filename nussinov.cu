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
        else cell_value = triInd_safe(i, j-1, len, DP);
#else
        cell_value = triInd_safe(i, j-1, len, DP);
#endif 

        // iterates through possible pairs for a cell
        for(int t = i; t< j-MIN_LOOP_LENGTH; t++){
          // Check paired scores (if pairing exists)
          if (pair_check(seq, t, j)) {
            int pairing1 = 0;
            int pairing2 = 0;

            if(i < (t-1)-MIN_LOOP_LENGTH) 
              pairing1 = triInd_safe(i, t-1, len, DP);
            if(t+1 < (j-1)-MIN_LOOP_LENGTH) 
              pairing2 = triInd_safe(t+1, j-1, len, DP);

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

void nussinov_gpu_wrap(uint8_t* seqs, uint32_t* seq_offsets, uint32_t* seq_lengths, uint32_t* dp_offsets, 
                        uint32_t N, uint32_t total_bytes, uint32_t total_dp_cells, int* batched_DP) {

  uint8_t* d_seqs;
  uint32_t *d_seq_offsets, *d_seq_lengths, *d_dp_offsets;
  int *d_batched_DP;

  Timer timerCudaMems;
  Timer timerCudaExec;
  Timer timerCudaFree;
  long time1, time2, time3;

  timerCudaMems.Start();
  CUDA_CHECK(cudaMalloc(&d_seqs, total_bytes * sizeof(uint8_t)));
  CUDA_CHECK(cudaMalloc(&d_seq_offsets, N * sizeof(uint32_t)));
  CUDA_CHECK(cudaMalloc(&d_seq_lengths, N * sizeof(uint32_t)));
  CUDA_CHECK(cudaMalloc(&d_dp_offsets, N * sizeof(uint32_t)));
  CUDA_CHECK(cudaMalloc(&d_batched_DP, total_dp_cells * sizeof(int)));

  CUDA_CHECK(cudaMemcpy(d_seqs, seqs, total_bytes * sizeof(uint8_t), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_seq_offsets, seq_offsets, N * sizeof(uint32_t), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_seq_lengths, seq_lengths, N * sizeof(uint32_t), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_dp_offsets, dp_offsets, N * sizeof(uint32_t), cudaMemcpyHostToDevice));
  
  // Ideally not necessary because all indexed values go on previous values or are 0
  //CUDA_CHECK(cudaMemset(d_batched_DP, 0, total_dp_cells * sizeof(int))); //calloc
  time1 = timerCudaMems.Stop();

  int num_blocks = N;

  timerCudaExec.Start();
  nussinov_gpu<<<num_blocks, BLOCK_SIZE>>>(d_seqs, d_seq_offsets, d_seq_lengths, d_dp_offsets, d_batched_DP, N);
  time2 = timerCudaExec.Stop();

  timerCudaFree.Start();
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());
  CUDA_CHECK(cudaMemcpy(batched_DP, d_batched_DP, total_dp_cells * sizeof(int), cudaMemcpyDeviceToHost));


  /*
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
  */
  CUDA_CHECK(cudaFree(d_seqs));
  CUDA_CHECK(cudaFree(d_seq_offsets));
  CUDA_CHECK(cudaFree(d_seq_lengths));
  CUDA_CHECK(cudaFree(d_dp_offsets));
  CUDA_CHECK(cudaFree(d_batched_DP));
  time3 = timerCudaFree.Stop();
  fprintf(stderr, "GPU exec in Memcpy: %ld Exec: %ld Free %ld\n\n", time1, time2, time3);
}
