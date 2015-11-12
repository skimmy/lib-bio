#include "common.h"

#include <iostream>

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

size_t prefixSuffixHammingDistance(const std::string& s1, const std::string& s2, size_t k) {
  //  std::cout << '\n' << s1 << '\t' << s2 << " (" << k << ")\t";
  size_t d = 0;
  size_t md = s2.length() - k;
  for (int i = 0; i < k; ++i, ++md) {
    //    std::cout << s1[i] << s2[md] << ' ';
    if (s2[i] != s1[md]) {
      d++;
    }
  }
  //  std::cout << "\t\t> " << d << " <\n";
  return d;
}
