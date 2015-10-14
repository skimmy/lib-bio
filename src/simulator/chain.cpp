#include "common.h"

#include <cstdint>
#include <iostream>

// half of this mxm matrix is not used in case we need some optimizations
uint64_t** M;

void initChainMatrix() {
  M = new uint64_t*[Options::opts.m];
  for (size_t i = 0; i <  Options::opts.m; ++i) {
    M[i] = new uint64_t[Options::opts.m];
    for (size_t j = 0; j < Options::opts.m; ++j) {
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
    for (size_t j = 0; j < Options::opts.m; ++j) {
      std::cout << M[i][j] << " ";
    }
    std::cout << '\n';
  }
}
