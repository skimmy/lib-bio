#ifndef DFT_H
#define DFT_H

#include <complex>
#include <vector>

#include "../core.h"

namespace dft {
  // functions to compute Discrete Fourier Transofmr
  typedef std::complex<double> tComp;
  std::vector<tComp> getUnityRoots(size_t n);
  tComp dftComponent(const std::vector<tComp>& in, size_t k, size_t n);

  // functions to encode bases into complex
  tComp baseToComplex(const char c);
  // encode a whole sequence into a vector of complex numbers
  std::vector<tComp> basesToComplexVector(const Sequence& s);

  /**
   * Change the map from base to complex numbers. It takes vector with 4
   * elements associated to the bases A,C,G and T respectively.
   */
  //void setEncodeVector(const std::vector<tComp>& newEnc);
}

#endif
