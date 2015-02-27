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

}
