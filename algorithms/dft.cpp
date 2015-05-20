#include "dft.hpp"

#include <complex>
#include <vector>
#include <cctype>


double const pi = 4.0 * std::real(std::atan(1.0));
double const sqrt2 = std::sqrt(2);


namespace dft {

  tComp encodedA(1,0);
  tComp encodedC(0,1);
  tComp encodedG(-1,0);
  tComp encodedT(0,-1);


  /**
   * Returns the n complex roots of unity
   */
  std::vector<tComp> getUnityRoots(size_t n) {
    tComp i(0,1);
    std::vector<tComp> omegas(n);
    for (size_t k = 0; k < n; ++k) {
      omegas[k] = std::polar(1.0,2*pi*k/n);
    }
    return omegas;
  }

  /**
   * Computes the k-th component of the order n DFT
   */
  tComp dftComponent(const std::vector<tComp>& in, size_t k, size_t n) {
    std::vector<tComp> omegas = getUnityRoots(n);
    tComp omega_k = std::pow(omegas[1],k);
    tComp result(0,0);
    for (size_t j = 0; j < n; ++j) {
      result += (omega_k *omegas[j]) * in[j];
    }
    return result;
  }

  /**
   * Transoform one base into
   */
  tComp baseToComplex(const char c) {
    switch (toupper(c)) {
    case 'A':
      return encodedA;
      break;
    case 'C':
      return encodedC;
      break;
    case 'G':
      return encodedG;
      break;
    case 'T':
      return encodedT;
      break;
    default:
      return tComp(0,0);
    }
  }

  std::vector<tComp> basesToComplexVector(const Sequence& s) {
    size_t n = s.getSequenceLength();
    std::vector<tComp> v(n);
    for (size_t k = 0; k < n; ++k) {
      v[k] = baseToComplex(s.getBaseAt(k));
    }
    return v;
  }

}
