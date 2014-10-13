#ifndef DNA_ALPHABET_2BITS_H
#define DNA_ALPHABET_2BITS_H

#include <map>

class DNAAlphabet2Bits {
public:
  static char intToChar(uint64_t i);
  static uint64_t charToInt(char c);

  static map<char, uint64_t> charToIntMap();
};

#endif
