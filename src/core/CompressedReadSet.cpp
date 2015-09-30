#include "CompressedReadSet.h"

#include <math.h>
#include <string.h>

#include <iostream>
#include <fstream>
using namespace std;

/***************** STATIC INITIALIZATIONS *******************/

int CompressedReadSet::BaseBitLength = 4;
int CompressedReadSet::DefaultReadSize = 48;


/********************** CONSTRUCTORS ************************/

CompressedReadSet::CompressedReadSet() {
  this->readCount = 0;
  this->readLength = CompressedReadSet::DefaultReadSize;
  this->elSize = CompressedReadSet::BaseBitLength;
  this->n = 0;
  init();
}

CompressedReadSet::CompressedReadSet(size_t readSize) {
 this->readCount = 0;
  this->readLength = readSize;
  this->elSize = CompressedReadSet::BaseBitLength;
  this->n = 0;
  init();
}

CompressedReadSet::CompressedReadSet(const Read* reads, size_t n) {
  this->readCount = n;
  if (this->readCount > 0) {
    this->elSize = CompressedReadSet::BaseBitLength;
    size_t m = (size_t) ceil( 8.0 / ((double)elSize) );
    this->readLength = (size_t) m * ceil((double)reads[0].length() / (double)m);
    this->n = this->readCount * this->readLength;
    init();
    memset(this->seq, 0, realSize);
    size_t k = 0;
    for (size_t i = 0; i < this->readCount; i++) {
      for (size_t j = 0; j < this->readLength; j++) {
	this->setElementAt(k++, DNACompressedSymbol::IupacToNumber(reads[i].getBases()[j]));
      }
    }
  }
}

CompressedReadSet::CompressedReadSet(const vector<Read>& reads) {
  this->readCount = reads.size();
  if (this->readCount > 0) {
    this->elSize = CompressedReadSet::BaseBitLength;
    size_t m = (size_t) ceil( 8.0 / ((double)elSize) );
    this->readLength = (size_t) m * ceil((double)reads[0].length() / (double)m);
    this->n = this->readCount * this->readLength;
    init();
    memset(this->seq, 0, realSize);
    // here copy content of 'reads' in sequence
    size_t k = 0;
    for (size_t i = 0; i < this->readCount; i++) {
      for (size_t j = 0; j < this->readLength; j++) {
	this->setElementAt(k++, DNACompressedSymbol::IupacToNumber(reads[i].getBases()[j]));
      }
    }
  }
}

CompressedReadSet::CompressedReadSet(const CompressedReadSet& other) 
  : CompressedSequence(other) {
  this->readCount = other.readCount;
  this->readLength = other.readLength;
}

CompressedReadSet::~CompressedReadSet() {
  if (this->seq) {
    delete[] this->seq;
    this->seq = 0;
  }
}

/******************* GET AND SET METHODS ********************/

size_t CompressedReadSet::getReadCount() const {
  return this->readCount;
}

uint8_t* CompressedReadSet::getReads(size_t start, size_t count) const {
  if ( start >= this->readCount ) {
    return 0;
  }
  uint8_t* retSeq = 0;
  size_t readRealSize = (size_t)ceil( (double)(this->readLength) / (8.0 / (double)(CompressedReadSet::BaseBitLength)));
  size_t m = min(count, this->readCount - start );
  size_t N = count * readRealSize;
  retSeq = new uint8_t[N];
  memset(retSeq, 0, N * sizeof(uint8_t));
  memcpy(retSeq, &(this->seq[readRealSize*start]), m * readRealSize * sizeof(uint8_t));
  return retSeq;
}

/********************* MODIFY METHODS ***********************/

CompressedReadSet& CompressedReadSet::append(const Read& read) {
  // adjust sequence size and read count
  this->readCount++;
  this->n += this->readLength;
  this->resize(this->n);
  int k = readLength * ( readCount - 1 );
  for (size_t j = 0; j < this->readLength; ++j) {
    this->setElementAt(k++, DNACompressedSymbol::IupacToNumber(read.getBases()[j]));
  }
  return *this;
}

CompressedReadSet& CompressedReadSet::append(const vector<Read>& reads) {
  size_t oldReadCount = this->readCount;
  this->readCount+= reads.size();
  this->n = this->readLength * this->readCount;
  this->resize(this->n);
  
  int k = readLength * oldReadCount;
  for (size_t i = 0; i < reads.size(); i++ ) {
    for (size_t j = 0; j < this->readLength; ++j) {
      this->setElementAt(k++, DNACompressedSymbol::IupacToNumber(reads[i].getBases()[j]));
    }
  }
  return *this;
}

/*********************** I/O METHODS ************************/

void CompressedReadSet::writeToFile(const string& fileName) const {
  ofstream ofs(fileName,ofstream::binary | ofstream::out);
  // write the real size of the array and size of an element
  ofs.write((char*)&(this->realSize), sizeof(size_t));
  ofs.write((char*)&(this->elSize), sizeof(size_t));
  ofs.write((char*)&(this->readCount), sizeof(size_t));
  ofs.write((char*)&(this->readLength), sizeof(size_t));
  // write the data
  ofs.write( (char*)this->seq, sizeof(uint8_t) * realSize );
}

CompressedReadSet& CompressedReadSet::loadFromFile(const string& fileName) {
  ifstream ifs(fileName, ofstream::binary | ofstream::in);
  // delete old sequence (if any)
  if (this->seq) {
    delete[] seq;
    seq = 0;
  }
  // load size information
  ifs.read((char*)&this->realSize, sizeof(size_t));
  ifs.read((char*)&this->elSize, sizeof(size_t));
  ifs.read((char*)&this->readCount, sizeof(size_t));
  ifs.read((char*)&this->readLength, sizeof(size_t));
  this->n = this->readCount * this->readLength;
  init();
  // load data
  ifs.read((char*)this->seq, this->realSize);
  return *this;
}

/************************************************************/
