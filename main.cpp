/* ECE213 Final Project
 * Authors: Francisco Gutierrez, Jaehyung Lim, Sean Lipps
 *
 * Code Nussinov Algorithm
 * Reference Python: https://github.com/cgoliver/Nussinov
 *
 * Date of Creation: 2/24/26
 * Date Last Modified: 2/25/26
 */



#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <bits/stdc++.h> // max function

#include "nussinov.cuh"


// cell index
struct cell_ind{
  int x;
  int y;
};

void show_DP(int* DP, int N){
  printf("Showing DP scores: \n");
  for(int i = 0; i<N; i++){
    for(int j = 0; j<N; j++){
      if(i < j - MIN_LOOP_LENGTH ) {
        printf("%3d ", DP[triInd(i, j, N)] );
      } else {
        printf("    ");
      }
    }
    printf("\n");
  }
}

void traceback(int i, int j, cell_ind* structure, int* DP, const uint8_t* seq, int* trace_len, int N) {
  if (j<=i){
    return;
  } else if ( DP[N*i+ j] == DP[N*i+ j-1] ){
    traceback(i, j-1, structure, DP, seq, trace_len, N);
  } else {
    for(int k = i; k < j-MIN_LOOP_LENGTH; k++) {
      if(pair_check(seq, k, j)){
        if (k-1<0) {
          if( DP[N*i+ j] ==  DP[N*(k+1)+ j-1] + 1) {
            structure[*trace_len].x = k;
            structure[*trace_len].y = j;
            (*trace_len)++;
            traceback(k+1, j-1, structure, DP, seq, trace_len, N);
            break;
          }
        } else if ( DP[N*i+j] == (
                      DP[N*(i)+k-1]
                    + DP[N*(k+1)+ j-1]
                    + 1
                  )) {
          // add the pair (j,k) to our list of pairs
          structure[*trace_len].x = k;
          structure[*trace_len].y = j;
          (*trace_len)++;
          // move the recursion to two substructures formed by this pairing
          traceback(i, k-1, structure, DP, seq, trace_len, N);
          traceback(k+1, j-1, structure, DP, seq, trace_len, N);
          break;
        }
      }
    }
  }
}

void write_structure(int N, cell_ind* structure, int* struct_len){
  char* dot_bracket = (char*)malloc(2*N+1);
  dot_bracket[N] = '\0';
  for(int i = 0; i<N; i++)
    dot_bracket[i]='.';

  for(int i = 0; i < *struct_len; i++){
    dot_bracket[std::min(structure[i].x, structure[i].y)] = '(';
    dot_bracket[std::max(structure[i].x, structure[i].y)] = ')';
  }

  bool put_comma = false;
  printf("('%s', [", dot_bracket);
  for(int i = 0; i < *struct_len; i++){
    if(put_comma) printf(", ");
    printf("(%d, %d)", structure[i].x, structure[i].y);
    put_comma = true;
  }
  printf("])\n");
  
  free(dot_bracket);
}

void nussinov_cpu(uint8_t* seq, int* DP, int N){
  int cell_value;
  int j;

  for(int k = MIN_LOOP_LENGTH+1; k < N; k++){
    for(int i = 0; i < (N-k); i++){
      j = i+k;

      if (i >= j- MIN_LOOP_LENGTH)
        cell_value = 0;
      else {
#if MIN_LOOP_LENGTH == 0
        // we do not want to go out of bounds
        if(j==0) cell_value = 0;
        else cell_value = DP[triInd(i, j-1, N)];
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

            cell_value = std::max(cell_value, pairing1 + pairing2 + 1);
          }
        }
      }

      DP[triInd(i, j, N)] = cell_value;
    }
  }
}

void nussinov(uint8_t* seq, int N){
  cell_ind *structure; // array tracing
  int* struct_len;
  int d_struct_len;
  struct_len = &d_struct_len; // may just want to ommit the pointer all together

  structure = (cell_ind *)malloc(2*N*sizeof(cell_ind)); 
  
  // annoying but needed for traceback
  int* DP_square = (int*)malloc( N*N*sizeof(int));

#ifdef CPU_TARGET

  int* DP = (int*)malloc( (((N-MIN_LOOP_LENGTH)*(N-MIN_LOOP_LENGTH-1)) /2 )*sizeof(int));

  nussinov_cpu(seq, DP, N);

#ifdef DEBUG
  show_DP(DP, N);
#endif

  // we allocate in bulk since we do not know the final traceback size
  // Using vectors may hurt us since these types do not exist and are 
  // not transferrable to GPUs
  *struct_len = 0;

  // Copy uptriangular matrix to real NxN Matrix
  for(int k = MIN_LOOP_LENGTH+1; k < N; k++){
    for (int i = 0; i < N-k; i++){
      int j = i+k;
      DP_square[N*i+j] = DP[triInd(i, j, N)];
    }
  }

  // uses square for traceback
  traceback(0, N-1, structure, DP_square, seq, struct_len, N);

  write_structure(N, structure, struct_len);

  printf("Running again on GPU\n");

#endif // CPU_TARGET

  nussinov_gpu_wrap(seq, DP_square, N);

  *struct_len = 0;

  traceback(0, N-1, structure, DP_square, seq, struct_len, N);

  write_structure(N, structure, struct_len);

#ifdef CPU_TARGET
  free(DP);
#endif
  free(DP_square);
  free(structure);
}

int main(int argc, char * const argv[]) {
  // nussinov sequence

  if (argc != 2) {
    printf("./nussinov <seq>\n");
    return 1;
  }

  int N = strlen(argv[1]);
  int num_bytes = (N+3)/4;
  uint8_t *seq = (uint8_t *)calloc(num_bytes, sizeof(uint8_t));
  if (!seq) return 1;

  // Pack the string into a 2-bit array
  for(int i = 0; i < N; i++) {
    uint8_t val = -1;
    switch(argv[1][i]) { // A is 00
        case 'A': case 'a': val = 0; break; // 00
        case 'C': case 'c': val = 1; break; // 01
        case 'G': case 'g': val = 2; break; // 10
        case 'T': case 't': val = 3; break; // 11
        case 'U': case 'u': val = 3; break; // 11
    }
    seq[i / 4] |= (val << ((i % 4) * 2)); // pack 2 bits into a byte
  }

#ifdef DEBUG
  printf("Length of sequence = %d\n", N);
#endif

  nussinov(seq, N);
}
