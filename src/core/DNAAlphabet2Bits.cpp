#include "../core.h"

map<char, uint64_t> initBasesMap() {
  map<char, uint64_t> basesMap;
  char basesUpper[] = "ACGT";
  char basesLower[] = "acgt";
  for (uint64_t i = 0; i < 4; ++i) {
    basesMap[basesUpper[i]] = i;
    basesMap[basesLower[i]] = i;
  }
  return basesMap;
}

// -----------------------------------------------------------------------------
//                               UTILITY FUNCTIONS
// -----------------------------------------------------------------------------
uint64_t fromChar(const char c) {
  static map<char, uint64_t> basesMap = initBasesMap();
  return basesMap[c] & 0x3;
}

char fromInt(const uint64_t i) {
  static char bases[] = "ACGT";
  return bases[(i & 0x3)];
}

// -----------------------------------------------------------------------------
//                               STATIC FUNCTIONS
// -----------------------------------------------------------------------------
char DNAAlphabet2Bits::intToChar(uint64_t i) {
  return fromInt(i); 
}

uint64_t DNAAlphabet2Bits::charToInt(char c) {
  return fromChar(c);
}

map<char, uint64_t> DNAAlphabet2Bits::charToIntMap() {
  return initBasesMap();
}
