#ifndef NUMERIC_KMER_H
#define NUMERIC_KMER_H

#include <core/Sequence.h>
#include <iostream>

namespace seq
{

// At the moment I don't feel this should extend the Sequence type,
// however to be on the safe side I will use same methods name as
// Sequence class
class NumericKMer {
private:
  uint64_t kmer;
  size_t k;

  // since it may be used many times conversion from string to int is
  // implemented in a separated (private) method. It also
  // automatically applies the mask.
  uint64_t sequenceToInt(const Sequence& s);

public:
  // CONSTRUCTORS
  NumericKMer(const Sequence& s);
  NumericKMer(size_t k);
  NumericKMer(uint64_t aKmer, size_t k);
  NumericKMer(const NumericKMer& aKmer);

  // CONVERSIONT OPERATORS
  operator uint64_t();

  // IOSTREAM
  friend std::ostream& operator<< (std::ostream& os, const NumericKMer& kmer);


  static const size_t KMax = 32;

  // STATIC UTILITY METHODS
  static uint64_t fromChars(const char* chars, size_t k);
};

}
#endif
