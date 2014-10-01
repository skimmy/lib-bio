#include "DNACompressedSymbol.h"

#include <map>
using namespace std;

map<char,int> initIupacMap() {
  map<char,int> iupacMap;
  char iupacsUp[] = "-ACMGRSVTWYHKDBN";
  char iupacsLo[] = ".acmgrsvtwyhkdbn";
  for (int i = 0; i < 16; i++) {
    iupacMap[iupacsUp[i]] = i;
    iupacMap[iupacsLo[i]] = i;
  }
  return iupacMap;
}

/********************** CONSTRUCTOR(S) **********************/

DNACompressedSymbol::DNACompressedSymbol(char c) {
  this->symbol = DNACompressedSymbol::IupacToNumber(c);
}

/************************ OPERATORS *************************/

DNACompressedSymbol& DNACompressedSymbol::operator = (const char c) {
  this->symbol = DNACompressedSymbol::IupacToNumber(c);
  return *this;
}

DNACompressedSymbol::operator char() const {
  return DNACompressedSymbol::NumberToIupac(this->symbol);
}

DNACompressedSymbol::operator uint8_t() const {
  return this->symbol;
}



/****************** STATIC UTILITY METHODS ******************/

char DNACompressedSymbol::NumberToIupac(uint8_t n) {
  static char iupacs[] = "-ACMGRSVTWYHKDBN";
  return iupacs[n];
}

uint8_t DNACompressedSymbol::IupacToNumber(char c) {
  static map<char,int> iupacMap = initIupacMap();
  return iupacMap[c];
}

/************************************************************/

