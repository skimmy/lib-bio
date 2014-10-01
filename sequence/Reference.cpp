#include "Reference.hpp"

#include "ColorAlphabet.hpp"

// ---------------------------------------------------------
//                 CONSTRUCTORS AND DESTRUCTOR
// ---------------------------------------------------------
Reference::Reference() {
  init(0);
}

Reference::Reference(const Reference& other) {
  init(other.length);
  memcpy(this->sequence, other.sequence, this->length);
}

Reference::Reference(const string& sequence) {
  init(sequence.length());
  memcpy(this->sequence, sequence.c_str(), this->length);
}

Reference::Reference(const char* sequence, size_t n) {
  init(n);
  memcpy(this->sequence, sequence, this->length);
}

Reference::~Reference() {
  if(this->sequence) {
    delete[] this->sequence;
  }
}

// ---------------------------------------------------------
//                    GET AND SET METHODS
// ---------------------------------------------------------
list<KMer> Reference::getKMerList(size_t k) {
  list<KMer> kmers;
  size_t m = this->length - k + 1;
  for (size_t i = 0; i < m; ++i) {
    kmers.push_back(KMer( &(this->sequence[i]) ,k));
  }
  return kmers;
}

// ---------------------------------------------------------
//                     CONVERSION METHODS
// ---------------------------------------------------------
Reference& Reference::toColors(char primer) {
  // TODO: insert check about sequence coding (i.e. isBaseSpace())
  string bases(this->sequence);
  string colors("");
  ColorAlphabet::basesToColors(bases, colors, primer);
  memcpy(this->sequence, colors.c_str(), this->length * sizeof(char));
  return *this;
}

Reference& Reference::toBases(char primer) {
  // TODO: insert check about sequence coding (i.e. isColorSpace())
  string colors(this->sequence);
  string bases("");
  ColorAlphabet::colorsToBases(colors, bases, primer);
  memcpy(this->sequence, bases.c_str(), this->length * sizeof(char));
  return *this;
}

// ---------------------------------------------------------
//                 'SEQUNCE' OVERRIDE METHODS
// ---------------------------------------------------------
const void* Reference::getSequence() const {
  return this->sequence;
}

size_t Reference::getSequenceLength() const {
  return (this->length * sizeof(char));
}
size_t Reference::getElementSize() const {
  return sizeof(char);
}
size_t Reference::getByteCount() const {
  return (length * sizeof(char));
}
char Reference::getBaseAt(size_t i) const {
  return this->sequence[i];
}

// ---------------------------------------------------------
//                      UTILITY METHODS
// ---------------------------------------------------------
void Reference::init(size_t n) {
  this-> length = n;
  this->sequence = new char[n];
}
