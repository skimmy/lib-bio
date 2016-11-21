#include "common.hpp"

#include <iostream>
#include <memory>

void
editDistanceEstimations(size_t n_min, size_t n_max, size_t n_step, size_t k_max) {
  
  std::cout << std::endl;
  for (size_t n = n_min; n <= n_max; n += n_step) {
    std::string s1(n,'N');
    std::string s2(n,'N');
    double AED = 0;
    for (size_t k = 1; k <= k_max; ++k) {
      generateIIDString(s1);
      generateIIDString(s2);
      AED += editDistance(s1,s2);
    }
    std::cout << n << "\t" << ( AED / k_max) << std::endl;
  }
  std::cout << std::endl;
}


std::unique_ptr<size_t[]>
editDistSamples(size_t n, size_t k_samples) {
  std::unique_ptr<size_t[]> v(new size_t[k_samples]);
  std::string s1(n,'N');
  std::string s2(n,'N');
  size_t* v0 = new size_t[n];
  size_t* v1 = new size_t[n];
  for (int k = 0; k < k_samples; ++k) {
    generateIIDString(s1);
    generateIIDString(s2);
    v[k] = editDistanceLinSpace(s1,s2,v0,v1);
  }
  delete[] v0;
  delete[] v1;
  return v;
}
