#include "hamming.hpp"

namespace bio {

  size_t hammingDistance(const Sequence& s1, const Sequence& s2) {
    size_t d = 0;
    size_t n = s1.getSequenceLength();
    for (size_t i = 0; i < n; ++i) {
      if (s1.getBaseAt(i) != s2.getBaseAt(i)) {
	d++;
      }
    }
    return d;
  }

}
