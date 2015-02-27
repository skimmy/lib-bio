#include "../algorithms.h"

namespace dstats
{

uint64_t computeD2(uint64_t X[], uint64_t Y[], int k) {
  size_t K = 1 << (2 * k);
  uint64_t D2 = 0;
  for (size_t i = 0; i < K; ++i) {
    D2 += (X[i] * Y[i]);
  }
  return D2;
}

uint64_t D2(const Sequence& x, const Sequence& y, size_t k) {
  size_t K = 1 << (2 * k);
  uint64_t* X = new uint64_t[K];
  uint64_t* Y = new uint64_t[K];
  spectrumAsArray(x,k,X);
  spectrumAsArray(y,k,Y);
  uint64_t D2 = computeD2(X,Y,k);
  delete[] Y;
  delete[] X;
  return D2;
}

}
