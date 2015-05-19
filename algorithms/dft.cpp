#include "dft.hpp"

#include <complex>
#include <vector>


double const pi = 4.0 * std::real(std::atan(1.0));

namespace dft {

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

}
