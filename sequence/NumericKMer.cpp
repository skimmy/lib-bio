#include <iostream>

#include "../sequence.h"

// -----------------------------------------------------------------------------
//                                  CONSTRUCTORS
// -----------------------------------------------------------------------------
NumericKMer::NumericKMer(const Sequence& s) {
  this->k = (s.getSequenceLength() < KMax) ? s.getSequenceLength() : KMax;  
  this->kmer = sequenceToInt(s);
}

NumericKMer::NumericKMer(size_t k) {
  this->kmer = 0x0;
  this->k = k;
}

NumericKMer::NumericKMer(uint64_t aKmer, size_t k) {
  this->kmer = aKmer;
  this->k = k;
}

NumericKMer::NumericKMer(const NumericKMer& aKmer) {
  this->kmer = aKmer.kmer;
  this->k = aKmer.k;
}

// -----------------------------------------------------------------------------
//                             CONVERSION OPERATORS
// -----------------------------------------------------------------------------

NumericKMer::operator uint64_t() {
  return this->kmer;
}


uint64_t NumericKMer::sequenceToInt(const Sequence& s) {
  uint64_t kmer = 0;
  //  uint64_t currentMask = 0x3;
  //  std::cout << "c  -   (2*i) base  -   kmer   -   mask" << std::endl;
  for (int i = 0; i < k; i++) {    
    char c = s.getBaseAt(i);
    // TODO: It may be worth to add a (fast) check for 'c'
    uint64_t base = DNAAlphabet2Bits::charToInt(c);
    //std::cout << std::hex << c << " (" << 2*i << ") " << base << " " << kmer << " " << currentMask << std::endl;
    kmer |= (base & 0x03);
    kmer <<= 2;
    //currentMask <<= 2;

  }
  return kmer;
}
