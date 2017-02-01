#include "generator.hpp"
#include "common.hpp"

#include <iostream>
#include <random>

namespace lbio { namespace sim { namespace  generator {

struct RandGen
{
  std::mt19937 gen;
  RandGen(std::mt19937 g) {
    gen = g;
    basesIdx = new std::uniform_int_distribution<>(0,3);
  }
  
  std::uniform_int_distribution<>* basesIdx;

  int nextBaseIdx() { return (*basesIdx)(gen); }
  
};

RandGen* randGen;

void
initRandomGenerator() {
  srand(time(NULL));
  std::random_device rd;
  std::mt19937 gen(rd());  
  randGen = new RandGen(gen);
}

int randomBaseIndex() {
  return randGen->nextBaseIdx();
}

char baseComplement(char b) {
  switch(b) {
  case 'A': case 'a':
    return 'T';

  case 'C': case 'c':
    return 'G';

  case 'G': case 'g':
    return 'C';

  case 'T': case 't':
    return 'A';    
  }
  return 'N';
}

char randomMutation(char c) { 
  return bases[(revBases[(int)c] + ((rand() % 3) + 1) ) & 0x3];
}

void simpleIIDErrors(std::string& s, double pe) {
  if (pe >0) {
    size_t m = s.length();
    for (size_t i = 0; i < m; ++i) {
      double x = (double)rand() / RAND_MAX;
      if (x < pe) {
	s[i] = randomMutation(s[i]);      
      }
    }
  }
}

void generateIIDGenome(size_t N, char* S) {		
  for (size_t i = 0; i < N; ++i) {		
    S[i] = bases[randomBaseIndex()];		
  }		
}

void generateIIDString(std::string& s) {
  for (size_t i = 0; i < s.size(); ++i) {
    s[i] = bases[randomBaseIndex()];
  }
}

void generateConstantGenome(size_t N, char* S, char b) {
  for (size_t i = 0; i < N; ++i) {
    S[i] = b;
  }
}

Read randomRead(size_t m) {
  char * tmp = new char[m];
  generateIIDGenome(m,tmp);
  Read r(std::string(tmp), 0);
  delete[] tmp;
  return r;
}

void simulateReadAt(size_t j, size_t m, const char* S, char* r) {
  for (size_t l = 0; l < m; ++l) {
    r[l] = S[j+l];
  }
}

void generateOfflineReads(const std::string& s, std::priority_queue<Read>& reads, size_t m, size_t M, double pe) {
  size_t N = s.length();
  const char* pS = s.c_str();
  size_t barN = N - m + 1;
  char* r = new char[m + 1]; // +1 accomodates null terminating symbol
  for (size_t i = 0; i < M; ++i) {
    size_t j = rand() % barN;
    simulateReadAt(j,m,pS,r);
    r[m] = 0;
    Read read(std::string(r),j);
    simpleIIDErrors(read.r, pe);
    reads.push(read);
  }
  delete[] r;
}

/**
 * Generates a noisy read starting at position j of S for online simulations
 */
Read generateOnlineRead(char* S, size_t j, size_t m, double pe) {  
  char* r = new char[m+1];
  simulateReadAt(j, m, S, r);
  Read read(std::string(r),j);
  simpleIIDErrors(read.r, pe);
  delete[] r;
  return read;
}

/**
 * Change all bases of a string substituting the original with the complement
 * (A <-> T C <-> G)
 */
void complementBases(char* S, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    S[i] = baseComplement(S[i]);
  }
}

} } } // namespaces
