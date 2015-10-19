#include "common.h"

#include <cstdint>
#include <iostream>

// half of this mxm matrix is not used in case we need some optimizations
uint64_t** M;

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

void evaluateChainRelation(const Read& r1, const Read& r2, size_t s) {
  if (s < Options::opts.m) {
    size_t e = prefixSuffixHammingDistance(r1.r, r2.r, s);
    //    std::cout << "(" << s << ", " << e << ")\t\t" << r1.r << " -- " << r2.r << "\n";
    M[s][e]++;
  }
}
