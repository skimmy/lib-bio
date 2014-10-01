#include "KMer.hpp"

// ---------------------------------------------------------
//                CONSTRUCTOS AND DESTRUCTOR
// ---------------------------------------------------------

KMer::KMer(size_t k) 
 : sequence("")
{
  this->k = k;
}

KMer::KMer(const string& seq) {
  sequence = seq;
  this->k = seq.length();
}

KMer::KMer(const char* seq, size_t k) {
  this->sequence = string(seq,k);
  this->k = k;
}

KMer::KMer(const KMer& other) {
  this->k = other.k;
  this->sequence = other.sequence;
}

// ---------------------------------------------------------
//                         OPERATORS
// ---------------------------------------------------------

ostream& operator<< ( ostream& os, const KMer& kmer) {
  os << kmer.sequence;
  return os;
}

KMer* KMer::operator=(const KMer& other) {
  this->sequence = other.sequence;
  this->k = other.k;
  return this;
}

bool KMer::operator==(const KMer& other) {
  return ( ( this->k == other.k ) && ( this->sequence == other.sequence ) );
}

// ---------------------------------------------------------
//                'SEQUENCE' OVERRIDE METHODS
// ---------------------------------------------------------

const void* KMer::getSequence() const {
  return this->sequence.c_str();
}

size_t KMer::getSequenceLength() const {
  return k;
}

size_t KMer::getElementSize() const {
  return sizeof(char);
}

size_t KMer::getByteCount() const {
  return k * sizeof(char);
}

char KMer::getBaseAt(size_t i) const {
  return sequence[i];
}

// ---------------------------------------------------------
