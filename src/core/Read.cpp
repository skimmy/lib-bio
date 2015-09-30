#include "Read.hpp"

using namespace std;

/*********************** CONSTRUCTORS ***********************/

Read::Read() 
  : header(""), bases(""), qualities("")
{
}

/******************** SET AND GET METHODS ********************/

void Read::setHeader(const string& header) {
  this->header = header;
}

void Read::setBases(const std::string& bases) {
  this->bases = bases;
}

void Read::setQualities(const std::string& qualities) {
  this->qualities = qualities;
}

std::string Read::getHeader() const {
  return this->header;
}

std::string Read::getBases() const {
  return this->bases;
}

std::string Read::getQualities() const {
  return this->qualities;
}

size_t Read::length() const {
  return bases.length();
}

std::list<KMer> Read::getKMerList(size_t k) const {
  std::list<KMer> kmersList;
  int M = bases.length() - k + 1;
  for (int i = 0; i < M; i++) {
    kmersList.push_back(KMer(bases.substr(i,k)));
  }
  return kmersList;
}

/******************** MODIFY OPERATIONS *********************/

Read& Read::trim(size_t n) {
  // TODO: check that n <= length(bases)
  bases.resize(n);
  qualities.resize(n);
  return *this;
}

/******************* 'SEQUENCE' OVERRIDE ********************/

const void* Read::getSequence() const {
  return this->bases.c_str();
}

size_t Read::getSequenceLength() const {
  return bases.length();
}

size_t Read::getElementSize() const {
  return sizeof(char);
}

size_t Read::getByteCount() const {
  return bases.length() * sizeof(char);
}

char Read::getBaseAt(size_t i) const {
  return this->bases[i];
}

/************************************************************/

