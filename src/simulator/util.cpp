#include <include/common.hpp>

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

// some util global variables (lookup tables, constants, ...)
double* power4_lookup = NULL; // contains 4^{-(m-s)} for s=0,...,m

void initUtil(size_t m) {
  power4_lookup = new double[m+1];
  for (size_t s = 0; s <= m; ++s) {
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

using Partition = std::vector<lbio_size_t>;
using ListOfPartitions = std::vector<Partition>;

ListOfPartitions
recursive_int_partition(lbio_size_t n, lbio_size_t k) {
  if (n == 0) {
    return ListOfPartitions();
  }
  if (k == 1) {
    Partition n_part(1,n);
    return ListOfPartitions(1, n_part);
  }
  ListOfPartitions P, tmp ;
  for (lbio_size_t p = 1; p <= n; p++) {
    tmp = recursive_int_partition(n-p, k-1);
    for (auto it_ = tmp.begin(); it_ != tmp.end(); it_++) {
      Partition P_p(1,p);      
      P_p.insert(P_p.end(), (*it_).begin(), (*it_).end());
      P.push_back(P_p);
    }
  }
  return P;     
}

