#include "DNAAlphabet.hpp"

const char SYMBOLS[] = {'A', 'C', 'G', 'T'};
const char SYMBOLS_LOWER_CASE[] = {'a', 'c', 'g', 't'};

map<char,size_t> DNAAlphabet::initCharMap() {
  map<char,size_t> outMap;
  for (int i = 0; i < 4; ++i) {
    outMap[SYMBOLS[i]] = i;
    outMap[SYMBOLS_LOWER_CASE[i]] = i;
  }
  return outMap;
}

size_t DNAAlphabet::length()  { 
  return LENGTH; 
}

char DNAAlphabet::getChar(size_t i, bool upper) {
  return upper ? SYMBOLS[i] : SYMBOLS_LOWER_CASE[i];
}

size_t DNAAlphabet::getIndex(char c) {
  static map<char,size_t> charMap = DNAAlphabet::initCharMap();
  return charMap[c];
}
