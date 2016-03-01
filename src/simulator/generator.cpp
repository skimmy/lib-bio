#include "common.hpp"

#include <iostream>

char randomMutation(char c) { 
  return bases[(revBases[c] + ((rand() % 3) + 1) ) & 0x3];
}

void simpleIIDErrors(std::string& s, double pe) {
  size_t m = s.length();
  for (size_t i = 0; i < m; ++i) {
    double x = (double)rand() / RAND_MAX;
    if (x < pe) {
      s[i] = randomMutation(s[i]);      
    }
  }
}

void generateIIDGenome(size_t N, char* S) {		
  for (size_t i = 0; i < N; ++i) {		
    S[i] = bases[rand() & 0x3];		
  }		
}

void generateConstantGenome(size_t N, char* S, char b) {
  for (size_t i = 0; i < N; ++i) {
    S[i] = b;
  }
}

void simulateReadAt(size_t j, size_t m, const char* S, char* r) {
  for (size_t l = 0; l < m; ++l) {
    r[l] = S[j+l];
  }
}

void generateOfflineReads(const std::string& s, std::priority_queue<Read>& reads) {
  size_t N = s.length();
  const char* pS = s.c_str();
  size_t m = Options::opts.m;
  size_t M = Options::opts.M;
  size_t barN = N - m + 1;
  char* r = new char[m + 1]; // +1 accomodates null terminating symbol
  for (size_t i = 0; i < M; ++i) {
    size_t j = rand() % barN;
    simulateReadAt(j,m,pS,r);
    r[m] = 0;
    Read read(std::string(r),j);
    simpleIIDErrors(read.r, Options::opts.pe);
    reads.push(read);
  }
  delete[] r;
}

/**
 * Generates a noisy read starting at position j of S for online simulations
 */
Read generateOnlineRead(char* S, size_t j) {  
  char* r = new char[Options::opts.m+1];
  simulateReadAt(j, Options::opts.m, S, r);
  Read read(std::string(r),j);
  delete[] r;
  return read;
}
