

#define MIN_LOOP_LENGTH 4

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 1
#endif

#ifndef GRID_SIZE
#define GRID_SIZE 1
#endif

void nussinov_gpu_wrap(char* seq, int* DP, int N);
