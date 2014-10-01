#include "FastqRead.h"

#define MAX_BUFFER_SIZE 4096

/********************** CONSTRUCTOR(S) **********************/

FastqRead::FastqRead() 
  : Read()
{
  
}

FastqRead::FastqRead(const FastqRead& other) {
  bases = other.bases;
  qualities = other.qualities;
  header = other.qualities;
}

/********************* STREAM OPERATORS *********************/

istream& operator>>(istream& is, FastqRead& read) {
  char buffer[MAX_BUFFER_SIZE];
  // header (starts with '@')
  is.getline(buffer,MAX_BUFFER_SIZE);
  read.header = string(buffer);
  // bases
  is.getline(buffer,MAX_BUFFER_SIZE);
  read.bases = string(buffer);
  //  '+' for qualities
  is.getline(buffer,MAX_BUFFER_SIZE);
  // qualtiies
  is.getline(buffer,MAX_BUFFER_SIZE);
  read.qualities = string(buffer);
  return is;
}

ostream& operator<<(ostream& os, const FastqRead& read) {
  os << read.header << '\n';
  os << read.bases << '\n';
  os << "+\n";
  os << read.qualities << '\n';
  return os;
}

/************************************************************/
