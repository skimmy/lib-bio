#include <core/CompressedSequence.h>

#include <iostream>
#include <fstream>
using namespace std;

#include <math.h>
#include <string.h>

/**************** CONSTRUCTOR(S) DESTRUCTOR *****************/

CompressedSequence::CompressedSequence() {
  this->n = 0;
  this->elSize = 0;
  this->mask = 0;
  this->realSize = 0;
  this->seq = 0;
}

CompressedSequence::CompressedSequence(size_t n, size_t elSize) {
  this->n = n;
  this->elSize = elSize;
  init();
  memset(this->seq, 0, this->realSize);
}

CompressedSequence::CompressedSequence(const CompressedSequence& s) {
  this->n = s.n;
  this->elSize = s.elSize;
  this->mask = s.mask;
  this->realSize = s.realSize;
  this->seq = new uint8_t[realSize];
  memcpy( this->seq, s.seq, this->realSize * sizeof(uint8_t));
}

CompressedSequence::~CompressedSequence() {
  delete[] seq;
}

/********************** QUERY METHODS ***********************/

size_t CompressedSequence::length() const {
  return this->realSize;
}

size_t CompressedSequence::getElementCount() const {
  return this->n;
}

uint8_t CompressedSequence::getElementAt(size_t i) const {
  return (seq[(i >> 1)] >> elSize*(1 - (i & 0x01))) & mask;
}

void CompressedSequence::setElementAt(size_t i, uint8_t e) {
  seq[(i >> 1)] |=  (mask & e) << elSize*(1 - (i & 0x01));
}

const uint8_t* CompressedSequence::getRawSequence() const {
  return this->seq;
}

/********************* MODIFY METHODS ***********************/

CompressedSequence& CompressedSequence::append(const CompressedSequence &other) {
  // size_t new_n = (this->n + other.n);
  size_t oldRealSize = this->realSize;
  this->realSize += other.realSize;
  this->n = this->n + other.n;
  uint8_t* new_sequence = new uint8_t[realSize];
  memcpy( new_sequence, this->seq, oldRealSize * sizeof(uint8_t));
  memcpy( &new_sequence[oldRealSize], other.seq, other.realSize * sizeof(uint8_t));
  delete[] seq;
  this->seq = new_sequence;
  return *this;
}

/************************ OPERATORS *************************/

// uint8_t CompressedSequence::operator[] (const size_t i) const {
//   return getElementAt(i);
// }

/*********************** I/O METHODS ************************/

void CompressedSequence::writeToFile(const string& fileName) const {
  ofstream ofs(fileName,ofstream::binary | ofstream::out);
  // write the real size of the array and size of an element
  ofs.write((char*)&(this->realSize), sizeof(size_t));
  ofs.write((char*)&(this->elSize), sizeof(size_t));
  // write the data
  ofs.write( (char*)this->seq, sizeof(uint8_t) * realSize );
}

CompressedSequence& CompressedSequence::loadFromFile(const string& fileName) {
  if (this->seq) {
    delete[] seq;
  }
  ifstream ifs(fileName, ofstream::in | ofstream::binary);
  // read real size of the array and size of an element
  ifs.read((char*)&(this->realSize), sizeof(size_t));
  ifs.read((char*)&(this->elSize), sizeof(size_t));
  init();
  // read the data
  ifs.read((char*)(this->seq), sizeof(uint8_t)*(this->realSize));
  return *this;
}

/******************* 'SEQUENCE' OVERRIDE ********************/

const void* CompressedSequence::getSequence() const {
  return this->seq;
}

size_t CompressedSequence::getSequenceLength() const {
  return this->n;
}

size_t CompressedSequence::getElementSize() const {
  return this->elSize;
}

size_t CompressedSequence::getByteCount() const {
  return this->realSize;
}

char CompressedSequence::getBaseAt(size_t i) const {
  return (char)getElementAt(i);
}

/********************* UTILITY METHODS **********************/

void CompressedSequence::init() {
  this->mask = 0xFF  >> ( 8*sizeof(uint8_t) - elSize );
  this->realSize = (size_t)ceil((double)n / (8.0 / (double)elSize));
  this->seq = new uint8_t[this->realSize];
}

void CompressedSequence::resize(size_t newSize) {
  size_t oldRealSize = this->realSize;
  this->n = newSize;
  uint8_t* temp = this->seq;
  init();
  memset(this->seq, 0, this->realSize);
  memcpy(this->seq, temp, oldRealSize * sizeof(uint8_t));
  delete[] temp;
}

/************************************************************/
