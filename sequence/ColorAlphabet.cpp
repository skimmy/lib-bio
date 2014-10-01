#include "ColorAlphabet.hpp"
#include "DNAAlphabet.hpp"

const char SYMBOLS[] = { '0', '1', '2', '3' };
const int TRANSLATION_MATRIX[4][4] = {
  { 0, 1, 2, 3 },
  { 1, 0, 3, 2 },
  { 2, 3, 0, 1 },
  { 3, 2, 1, 0 }
};

map<char,size_t> ColorAlphabet::initCharMap() {
  map<char,size_t> outMap; 
  for (int i = 0; i < 4; ++i) {
    outMap[SYMBOLS[i]] = i;
  }
  return outMap;
}


size_t ColorAlphabet::length() {
  return LENGTH;
}

char ColorAlphabet::getChar(size_t i) {
  return SYMBOLS[i];
}

size_t ColorAlphabet::getIndex(char c) {
  static map<char,size_t> charMap = ColorAlphabet::initCharMap();
  return charMap[c];
}

void ColorAlphabet::colorsToBases(const string& colors, string& bases, char primer) {
  bases.clear();
  size_t N = colors.size();
  size_t i = DNAAlphabet::getIndex(primer);
  size_t j = 0;
  for (size_t k = 0; k < N; ++k) {
    j = ColorAlphabet::getIndex(colors[k]);
    int nextBaseIndex =  TRANSLATION_MATRIX[i][j];
    bases += DNAAlphabet::getChar(nextBaseIndex);
    i = nextBaseIndex;
  }
}

void ColorAlphabet::basesToColors(const string& bases, string& colors, char primer) {
  colors.clear();
  size_t N = bases.size();
  size_t i = DNAAlphabet::getIndex(primer);
  size_t j = 0;
  for (size_t k = 0; k < N; ++k) {
    j = DNAAlphabet::getIndex(bases[k]);
    int nextColorIndex = TRANSLATION_MATRIX[i][j];
    colors += ColorAlphabet::getChar(nextColorIndex);
    i = j;
  }
}
