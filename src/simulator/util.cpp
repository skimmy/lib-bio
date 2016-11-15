#include "common.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

// some util global variables (lookup tables, constants, ...)
double* power4_lookup = NULL; // contains 4^{-(m-s)} for s=0,...,m

void initUtil() {

  int m = Options::opts.m;
  power4_lookup = new double[m+1];
  for (int s = 0; s <= m; ++s) {
    power4_lookup[s] = pow(4, s-m);
  }
}

void clearUtil() {  
  delete[] power4_lookup;
  power4_lookup = NULL;
}

double** initDoubleMatrix(size_t n, size_t m) {
  double** matrix = new double*[n];
  for (size_t i = 0; i < n; ++i) {
    matrix[i] = new double[m];
    for (size_t j = 0; j < m; ++j) {
      matrix[i][j] = 0.0;
    }
  }
    return matrix;
}

void clearDoubleMatrix(double** matrix, size_t n, size_t m) {
  for (size_t i = 0; i < n; ++i) {
    delete[] matrix[i];
  }
  delete[] matrix;
}

double elementsSumDoubleMatrix(double** matrix, size_t n, size_t m) {
  double sum = 0.0;
  for (size_t i  = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      sum += matrix[i][j];
    }
  }
  return sum;
}

void printString(char* s, size_t n) {
  for (int i = 0; i < n; ++i) {
    std::cout << s[i];
  }
}

void printDoubleMatrix(double** M, size_t n, size_t m) {
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      std::cout << M[i][j] << '\t';
    }
    std::cout << '\n';
  }
}

size_t hammingDistance(const char* s1, const char* s2, size_t m) {
  size_t d = 0;
  for (int i = 0; i < m; ++i) {
    if (s1[i] != s2[i]) {
      ++d;
    }
  }
  return d;
}

size_t hammingDistance(const std::string& s1, const std::string& s2, size_t m) {
  return hammingDistance(s1.c_str(), s2.c_str(), m);
}

// Calculates Hamming distance between k suffix of s1 and k prefix of s2 
size_t prefixSuffixHammingDistance(const std::string& s1, const std::string& s2, size_t k) {
  size_t d = 0;
  size_t md = s2.length() - k;
  for (int i = 0; i < k; ++i, ++md) {
    if (s2[i] != s1[md]) {
      d++;
    }
  }
  return d;
}

size_t bestHammingOverlap(const std::string& s1, const std::string& s2) {
  size_t m = s1.length();
  for (size_t i = m; i >= 1; --i) {
    size_t d_i = prefixSuffixHammingDistance(s2, s1, i);
    if (d_i == 0) {
      return i;
    }
  }
  return 0;
}

// Edit distance functions
size_t
editDistance(const std::string& s1, const std::string& s2) {
  size_t n = s1.size();
  size_t m = s2.size();
  size_t** dpMatrix = new size_t*[n+1];
  for (size_t i = 0; i < n+1; ++i) {
    dpMatrix[i] = new size_t[m+1];
  }

  // initialization of first row and column
  for (size_t i = 0; i < n+1; ++i) {
    dpMatrix[i][0] = i;
  }
  for (size_t j = 0; j < m+1; ++j) {
    dpMatrix[0][j] = j;
  }

  for (int i = 1; i < n+1; ++i) {
    for(int j = 1; j < m+1; ++j) {
      size_t delta = (s1[i-1] == s2[j-1]) ? 0 : 1;
      dpMatrix[i][j] = MIN( MIN(dpMatrix[i-1][j]+1, dpMatrix[i][j-1]+1) , dpMatrix[i-1][j-1] + delta ) ;
    }
  }
  size_t dist = dpMatrix[n][m];
  for (int i = 0; i < n+1; ++i) {
    delete[] dpMatrix[i];   
  }
  delete[] dpMatrix;
  return dist;
}

// Most of these function are used for debug and dev purpose, but they can be
// useful in other situations (perhaps in the future development of simulator)

Read randomRead(size_t m) {
  char * tmp = new char[m];
  generateIIDGenome(m,tmp);
  Read r(std::string(tmp), 0);
  delete[] tmp;
}

// Error

void
fatal_error(const std::string &msg, int exit_code)
{
  std::cerr << "[ERROR]  " + msg << std::endl;
  exit(exit_code);
}
