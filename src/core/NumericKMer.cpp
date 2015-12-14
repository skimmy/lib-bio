#include <iostream>

#include "../core.h"

namespace seq
{

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
  for (size_t i = 0; i < k; ++i) {    
    kmer <<= 2;
    char c = s.getBaseAt(i);
    // TODO: It may be worth to add a (fast) check for 'c'
    uint64_t base = DNAAlphabet2Bits::charToInt(c);
    kmer |= (base & 0x03);   
  }
  return kmer;
}

// -----------------------------------------------------------------------------
//                            STATIC UTILITY METHODS
// -----------------------------------------------------------------------------
uint64_t NumericKMer::fromChars(const char* chars, size_t k) {
  uint64_t kmer = 0;
  for (size_t i = 0; i < k; ++i) {
    kmer <<= 2;
    char c = chars[i];
    uint64_t base = DNAAlphabet2Bits::charToInt(c);
    kmer |= (base & 0x03);
  }
  return kmer;
}

// -----------------------------------------------------------------------------
//                                     IOSTREAM
// -----------------------------------------------------------------------------

ostream& operator<< (ostream& os, const NumericKMer& kmer) {
  uint64_t copy = kmer.kmer;
  size_t k = kmer.k;
  std::string out("");
  for (size_t i = 0; i < k; ++i) {
    out = DNAAlphabet2Bits::intToChar(copy & 0x03) + out;
    copy >>= 2;
 }	 
 os << out;
 return os;
}

}
