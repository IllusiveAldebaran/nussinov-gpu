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

#define MIN_LOOP_LENGTH 4

// cell index
struct cell_ind{
  int x;
  int y;
};

void show_DP(int* DP, int N){
  printf("Showing DP scores: \n");
  for(int i = 0; i<N; i++){
    for(int j = 0; j<N; j++){
      printf("%3d ", *((DP+N*i)+j) );
    }
    printf("\n");
  }
}


inline bool pair_check(const uint8_t* seq, int i, int j) {
  // retrieve 2 bits from a byte with mask
  uint8_t nuc1 = (seq[i / 4] >> ((i % 4) * 2)) & 0x03; 
  uint8_t nuc2 = (seq[j / 4] >> ((j % 4) * 2)) & 0x03;
  return (nuc1 ^ nuc2) == 3;
}

int* initialize(int N) {
  // NxN matrix with scores of optimal pairings
  int* DP = (int*)malloc(N*N*sizeof(int));
  int j;
  for(int k =0; k< MIN_LOOP_LENGTH; k++) {
    for(int i = 0; i < N-k; i++){
      j = i + k;
      *((DP+i*N)+j) = 0; //INT_MIN;
    }

  }

  return DP;
}

int opt(int i, int j, const uint8_t* seq){
  if (i>= j-MIN_LOOP_LENGTH){
    return 0;
  }

  int unpaired = opt(i, j-1, seq);

  // TODO: This is a reduction problem
  int paired = 0; // (maximum)
  for(int t = i; t< j-MIN_LOOP_LENGTH;t++){
    if (pair_check(seq, t, j)) {
      // TODO: This recursiveness HAS to be incredible inefficient...
      paired = std::max(paired, 1+opt(i, t-1, seq) + opt(t+1, j-1, seq));
    }
  }

  return std::max(unpaired, paired);
}

void traceback(int i, int j, cell_ind* structure, int* DP, const uint8_t* seq, int* trace_len, int N) {
  if (j<=i){
    return;
  } else if ( *((DP+N*i)+j) == *((DP+N*i)+(j-1)) ){
    traceback(i, j-1, structure, DP, seq, trace_len, N);
  } else {
    for(int k = i; k < j-MIN_LOOP_LENGTH; k++) {
      if(pair_check(seq, k, j)){
        if (k-1<0) {
          if( *((DP+N*i)+j) ==  *((DP+N*(k+1))+(j-1)) + 1) {
            structure[*trace_len].x = k;
            structure[*trace_len].y = j;
            (*trace_len)++;
            traceback(k+1, j-1, structure, DP, seq, trace_len, N);
            break;
          }
        } else if ( *((DP+N*i)+j) == (
                      *((DP+N*i)+(k-1))
                    + *((DP+N*(k+1))+(j-1))
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

void write_structure(int N, cell_ind* structure,int* struct_len){
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


void nussinov(uint8_t* seq, int N){
  cell_ind *structure; // array tracing
  int* struct_len;
  int d_struct_len;
  struct_len = &d_struct_len; // may just want to ommit the pointer all together

  int* DP = initialize(N);

  {
    int j;
    for(int k = MIN_LOOP_LENGTH; k < N; k++){
      for(int i = 0; i < (N-k); i++){
        j = i+k;
        *((DP+i*N)+j) = opt(i, j, seq);
      }
    }
  }

  // Copy values to lower triangle to avoid null references
  for(int i = 0; i < N; i++){
    for(int j = 0; j < i; j++){
      *((DP+N*i)+j) = *((DP+N*j)+i);
    }
  }


#ifdef DEBUG
  show_DP(DP, N);
#endif

  // we allocate in bulk since we do not know the final traceback size
  // Using vectors may hurt us since these types do not exist and are 
  // not transferrable to GPUs
  structure = (cell_ind *)malloc(2*N*sizeof(cell_ind)); 
  *struct_len = 0;

  traceback(0, N-1, structure, DP, seq, struct_len, N);

  write_structure(N, structure, struct_len);

  free(seq);
  free(DP);
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
    uint8_t val = 0;
    switch(argv[1][i]) { // A is 00
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
