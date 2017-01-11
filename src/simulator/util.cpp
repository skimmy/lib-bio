#include "common.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

// some util global variables (lookup tables, constants, ...)
double* power4_lookup = NULL; // contains 4^{-(m-s)} for s=0,...,m

char bases[] = {'A', 'C', 'G', 'T'};
char revBases[128];

void initUtil() {
  revBases['A'] = revBases['a'] = 0;
  revBases['C'] = revBases['c'] = 1;
  revBases['G'] = revBases['g'] = 2;
  revBases['T'] = revBases['t'] = 3;
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


double elementsSumDoubleMatrix(double** matrix, size_t n, size_t m) {
  double sum = 0.0;
  for (size_t i  = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      sum += matrix[i][j];
    }
  }
  return sum;
}

// FUNCTIONS FOR 2-BITS ENCODING CONVERSION
uint64_t string2Encode(const std::string&s) {
  size_t n = std::max<int>(32, s.size());
  uint64_t mask = 0;
  uint64_t enc = 0;
  for (size_t i = 0; i < n; ++i) {
    uint64_t b = revBases[(size_t)s[i]];
    enc = enc | ( (b & 0x3) << 2*i );
    mask = mask | (0x3 << 2*i);
  }
  return (enc & mask);
}

std::string encoding2String(uint64_t e, size_t n) {
  std::string s;
  for (size_t i = 0; i < n; ++i) {
    s += bases[e & 0x03];
    e = e >> 2;
  }
  return s;
}

size_t hammingDistance(const char* s1, const char* s2, size_t m) {
  size_t d = 0;
  for (size_t i = 0; i < m; ++i) {
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
  for (size_t i = 0; i < k; ++i, ++md) {
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


void printVec(size_t* v, size_t n) {
  for (size_t i =0; i <n;i++){
    std::cout << v[i] << ' ';
  }
  std::cout << std::endl;
}

// Most of these function are used for debug and dev purpose, but they can be
// useful in other situations (perhaps in the future development of simulator)

Read randomRead(size_t m) {
  char * tmp = new char[m];
  generateIIDGenome(m,tmp);
  Read r(std::string(tmp), 0);
  delete[] tmp;
  return r;
}

// Error and warnings


void
fatal_error(const std::string &msg, int exit_code)
{
  logError(msg);
  exit(exit_code);
}
