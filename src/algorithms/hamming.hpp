#ifndef HAMMING_H
#define HAMMING_H

#include <core/Sequence.h>

namespace bio {
  /**
   * Compute Hamming distance between given sequences. Note that sequences must
   * have same length or unexpected behavior will occur.
   */
  size_t hammingDistance(const Sequence& s1, const Sequence& s2);
}

#endif
