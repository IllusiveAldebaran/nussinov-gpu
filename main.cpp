/* ECE213 Final Project
 * Authors: Francisco Gutierrez, Jaehyung Lim, Sean Lipps
 *
 * Code Nussinov Algorithm
 * Reference Python: https://github.com/cgoliver/Nussinov
 *
 * Date of Creation: 2/24/26
 * Date Last Modified: 2/24/26
 */



#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <bits/stdc++.h> // max function

#include "nussinov.cuh"

// TODO Encode 
// ACGT: 
// A: 00
// C: 01
// G: 10
// T: 11 (U)

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

constexpr bool pair_check(char nuc1, char nuc2) {
  bool check = false;

  // TODO: encode and compare via bool
  if(nuc1 == 'A' && nuc2 == 'U') check = true;
  if(nuc1 == 'U' && nuc2 == 'A') check = true;
  if(nuc1 == 'C' && nuc2 == 'G') check = true;
  if(nuc1 == 'G' && nuc2 == 'C') check = true;

  return check;
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

int opt(int i, int j, char* seq){
  if (i>= j-MIN_LOOP_LENGTH){
    return 0;
  }

  int unpaired = opt(i, j-1, seq);

  // TODO: This is a reduction problem
  int paired = 0; // (maximum)
  for(int t = i; t< j-MIN_LOOP_LENGTH;t++){
    if (pair_check(seq[t], seq[j])) {
      // TODO: This recursiveness HAS to be incredible inefficient...
      paired = std::max(paired, 1+opt(i, t-1, seq) + opt(t+1, j-1, seq));
    }
  }

  return std::max(unpaired, paired);
}

void traceback(int i, int j, cell_ind* structure, int* DP, char* seq, int* trace_len, int N) {
  if (j<=i){
    return;
  } else if ( *((DP+N*i)+j) == *((DP+N*i)+(j-1)) ){
    traceback(i, j-1, structure, DP, seq, trace_len, N);
  } else {
    for(int k = i; k < j-MIN_LOOP_LENGTH; k++) {
      if(pair_check(seq[k], seq[j])){
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

void write_structure(char* seq, int N, cell_ind* structure,int* struct_len){
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


void nussinov(char* seq, int N){
  char *dot_bracket;
  cell_ind *structure; // array tracing
  int* struct_len;
  int d_struct_len;
  struct_len = &d_struct_len; // may just want to ommit the pointer all together

  int* DP = initialize(N);

  //{
  //  int j;
  //  for(int k = MIN_LOOP_LENGTH; k < N; k++){
  //    for(int i = 0; i < (N-k); i++){
  //      j = i+k;
  //      *((DP+i*N)+j) = opt(i, j, seq);
  //    }
  //  }
  //}

  //// Copy values to lower triangle to avoid null references
  //for(int i = 0; i < N; i++){
  //  for(int j = 0; j < i; j++){
  //    *((DP+N*i)+j) = *((DP+N*j)+i);
  //  }
  //}


#ifdef DEBUG
  show_DP(DP, N);
#endif

  // we allocate in bulk since we do not know the final traceback size
  // Using vectors may hurt us since these types do not exist and are 
  // not transferrable to GPUs
  structure = (cell_ind *)malloc(2*N*sizeof(cell_ind)); 
  //*struct_len = 0;

  //traceback(0, N-1, structure, DP, seq, struct_len, N);

  //write_structure(seq, N, structure, struct_len);


  //printf("Running again on GPU\n");

  nussinov_gpu_wrap(seq, DP, N);

  // Copy values to lower triangle to avoid null references
  for(int i = 0; i < N; i++){
    for(int j = 0; j < i; j++){
      *((DP+N*i)+j) = *((DP+N*j)+i);
    }
  }

  *struct_len = 0;

  traceback(0, N-1, structure, DP, seq, struct_len, N);

  write_structure(seq, N, structure, struct_len);

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

  char *seq;
  int N = 0;

  N = strlen(argv[1]);
#ifdef DEBUG
  printf("Length of sequence = %d\n", N);
#endif

  // TODO:
  // check for vulnerability
  // clean up the input
  if ((seq = (char *)malloc(N + 1)) != NULL) {
    bzero(seq, N + 1); 
    strncpy(seq, argv[1], N);
#ifdef DEBUG
    printf("argv[1] = %s\n", seq);
#endif
  }

  nussinov(seq, N);
}
