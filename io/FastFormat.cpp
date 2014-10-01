#include "FastFormat.h"

#include <iostream>
#include <fstream>

const int MAX_BUFFER_SIZE = 4096;

// ---------------------------------------------------------
//                       CONSTRUCTORS
// ---------------------------------------------------------
FastFormat::FastFormat()
  : Format("FAST"), sequence(""), header("")
{
}

// ---------------------------------------------------------
//                     CONVERSION METHODS
// ---------------------------------------------------------
Reference FastFormat::toReference() const {
  return (Reference)(*this);
}

// ---------------------------------------------------------
//                         OPERATORS
// ---------------------------------------------------------
FastFormat::operator Reference() const {
  return Reference(this->sequence);
}

// ---------------------------------------------------------
//                 'FORMAT' METHODS OVERRIDE
// ---------------------------------------------------------
string FastFormat::loadFromFile(const string &fileName) {
  ifstream ifs(fileName.c_str(), ifstream::in);
  char buffer[MAX_BUFFER_SIZE];
  while(!ifs.eof()) {
    ifs.getline(buffer,MAX_BUFFER_SIZE);
    if ( buffer[0] == '>' ) {
      header += string(buffer) + '\n';
      continue;
    }
    sequence += string(buffer);
  }
  return sequence;
}

string FastFormat::getSequence() const {
  return this->sequence;
}

string FastFormat::getHeader() const {
  return this->header;
}

// ---------------------------------------------------------
