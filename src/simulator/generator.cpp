#include "common.h"

void generateIIDGenome(size_t N, char* S) {		
  for (size_t i = 0; i < N; ++i) {		
    S[i] = bases[rand() & 0x3];		
  }		
}		
