#include "common.h"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <map>

// matrix of false positive
double ** fpMatrix;

// half of this mxm matrix is not used in case we need some optimizations
uint64_t** M;
// this map keeps track of the distance between non overlapping reads
std::map<size_t, size_t> D;

void initFalsePositiveMatrix() {

  double p_no_olap = (double)(Options::opts.N - 2 * Options::opts.m)
    / (double)pow(4,Options::opts.m);
  
  
  fpMatrix = new double*[Options::opts.m];
  for (size_t i = 0; i < Options::opts.m; ++i) {
    fpMatrix[i] = new double[Options::opts.m + 1];
    for (size_t j = 0; j < Options::opts.m + 1; ++j) {

      double p_olap = 0.0;
      double p_hat_s = 0.0;
      // NOTE: For now we use the equation for
      //   P[r2 = B | r1 = A, J = j, D = \Delta]
      // however this should (probably) be extended to all values of r1 = A

      fpMatrix[i][j] = (p_no_olap + p_olap + p_hat_s) / (double)Options::opts.N;
    }
  }
}

void clearFalsePositiveMatrix() {
   for (size_t i = 0; i <  Options::opts.m; ++i) {
    delete[] fpMatrix[i];
  }
  delete[] fpMatrix;
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
    //    std::cout << "(" << s << ", " << e << ")\t\t" << r1.r << " -- " << r2.r << "\n";
    M[s][e]++;
  } 
}

void addNonOverlapRecord(size_t d) {
  D[d]++;
}
