#include <include/common.hpp>
#include <include/chain.hpp>


#include <include/options.hpp>

#include <include/generator.hpp>
#include <include/util.hpp>


#include <cstdint>
#include <cmath>
#include <iostream>
#include <map>

using namespace lbio::sim::generator;

// matrix of false positive
double ** fpMatrix;
// half of this mxm matrix is not used in case we need some optimizations
uint64_t** M;
// this map keeps track of the distance between non overlapping reads
std::map<size_t, size_t> D;

void computeExpectedFalsePositiveMatrix(double** m) {
  for (size_t i = 0; i < Options::opts.m; ++i) {
    for (size_t j = 0; j < Options::opts.m + 1; ++j) {
      m[i][j] = fpMatrix[i][j] * (double)M[i][j];
      //      std::cout << fpMatrix[i][j] << ' ';
    }
  }
}

void printFalsePositiveMatrix() {
    std::cout << '\n';
  for (size_t i = 0; i <  Options::opts.m; ++i) {   
    for (size_t j = 0; j < Options::opts.m + 1; ++j) {
      std::cout << fpMatrix[i][j] << "\t";
    }
    std::cout << '\n';
  }
}


void initChainMatrix() {
  M = new uint64_t*[Options::opts.m];
  for (size_t i = 0; i <  Options::opts.m; ++i) {
    M[i] = new uint64_t[Options::opts.m + 1];
    for (size_t j = 0; j < Options::opts.m+1; ++j) {
      M[i][j] = 0;
    }
  }
}

void clearChainMatrix() {   
  for (size_t i = 0; i <  Options::opts.m; ++i) {
    delete[] M[i];
  }
  delete[] M;
}

void printChainMatrix() {
  std::cout << '\n';
  for (size_t i = 0; i <  Options::opts.m; ++i) {   
    for (size_t j = 0; j < Options::opts.m + 1; ++j) {
      std::cout << M[i][j] << "\t";
    }
    std::cout << '\n';
  }
}

void printNonOverlapDistribution() {
  for (std::pair<size_t,size_t> d : D) {
    std::cout << (int)d.first << "\t" << d.second << std::endl;
  }
}



void evaluateChainRelation(const Read& r1, const Read& r2, size_t s) {
  if (s < Options::opts.m) {
    size_t e = prefixSuffixHammingDistance(r1.r, r2.r, s);
    M[s][e]++;
  } 
}

void addNonOverlapRecord(size_t d) {
  D[d]++;
}
