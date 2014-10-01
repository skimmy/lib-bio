#include "FastqFormat.h"

#include <iostream>
#include <fstream>
using namespace std;

const int MAX_BUFFER_SIZE = 4096;

/******************** (DE)CONSTRUCTOR(S) ********************/

FastqFormat::FastqFormat() 
  : Format("FASTQ"), sequence(""), header(""), inFile()
{
  nextRead = 0;
}

FastqFormat::~FastqFormat() {
  if (nextRead) {
    delete nextRead;
  }
}

/*********************** LOAD METHODS ***********************/

string FastqFormat::loadFromFile(const string &fileName) {
  ifstream ifs(fileName.c_str(), ifstream::in);
  char buffer[MAX_BUFFER_SIZE];
  while(!ifs.eof()) {
    ifs.getline(buffer,MAX_BUFFER_SIZE);
    sequence += string(buffer) + '\n';
  }
  return sequence;
}

bool FastqFormat::openFile(const string& fileName) {
  inFile.open(fileName);
  return inFile.good();
}

bool FastqFormat::hasNextRead() {
  return (nextRead != 0 && nextRead->getBases().length() > 0);
}

FastqRead FastqFormat::getNextRead() {
  FastqRead toReturn = (nextRead) ? FastqRead(*nextRead) : FastqRead();
  loadNextRead();  
  return toReturn;
}



/*********************** GET METHODS ************************/

string FastqFormat::getSequence() const {
  return this->sequence;
}
string FastqFormat::getHeader() const {
  return this->header;
}

/******************** UTILITY  METHODS **********************/

//readsFile >> read;


void FastqFormat::loadNextRead() {
  if (nextRead) {
    delete nextRead;
  }
  nextRead = new FastqRead();
  if (inFile.is_open()) {
    inFile >> *nextRead;
  }
}

/************************************************************/
