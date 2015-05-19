#ifndef DFT_H
#define DFT_H

#include <complex>
#include <vector>

namespace dft {
  typedef std::complex<double> tComp;
  std::vector<tComp> getUnityRoots(size_t n);
  tComp dftComponent(const std::vector<tComp>& in, size_t k, size_t n);
}

#endif
