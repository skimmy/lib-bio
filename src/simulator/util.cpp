#include "common.h"

#include <iostream>

void printString(char* s, size_t n) {
  for (int i = 0; i < n; ++i) {
    std::cout << s[i];
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

size_t prefixSuffixHammingDistance(const std::string& s1, const std::string& s2, size_t k) {
  //  std::cout << '\n' << s1 << '\t' << s2 << " (" << k << ")\t";
  size_t d = 0;
  size_t md = s2.length() - k;
  for (int i = 0; i < k; ++i, ++md) {
    //std::cout << s1[i] << s2[md] << ' ';
    if (s1[i] != s2[md]) {
      d++;
    }
  }
  //  std::cout << '\n';
  return d;
}
