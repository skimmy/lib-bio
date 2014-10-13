#include <vector>
#include <map>

#include <iostream>

#include "../alignment.h"
#include "../sequence.h"

// char basesArray[] = "ACGT";
// std::map< char, uint64_t > constructBasesMap() {
//   std::map< char, uint64_t > basesMap;
//   basesMap['A'] = 0;
//   basesMap['C'] = 1;
//   basesMap['G'] = 2;
//   basesMap['T'] = 3;
//   return basesMap;
// }

// string codeToBases(uint64_t code, size_t k) {
//   string s = "";
//   for (size_t i = 0; i < k; ++i) {
//     s = basesArray[code & 0x3] + s;
//     code = code >> 2;
//   }
//   return s;
// }

std::vector<uint64_t>* getKmersFrequency(const Sequence& seq, size_t k) {
  std::ios::fmtflags f( std::cout.flags() );
  // TODO: check for k < 0 and k > N
  size_t N = seq.getSequenceLength();
  size_t M = (1 << 2*k);
  std::map< char, uint64_t > basesMap = DNAAlphabet2Bits::charToIntMap();
  // init the frequencies vector
  std::vector<uint64_t>* freq = new std::vector<uint64_t>(M);
  for (size_t l = 0; l < M; ++l) {
    (*freq)[l] = 0;
  }
  uint64_t mask = M - 1;

  // construct the code for the first kmer
  uint64_t index = basesMap[seq.getBaseAt(0)];
  for (size_t j = 1; j < k; ++j) {    
    // shift and insert the next code
    index = ((index << 2) & mask)  | (basesMap[seq.getBaseAt(j)] & 0x3);
  }

  ((*freq)[index])++;
  // scan for the next 
  for (size_t i = k; i < N; ++i) {
    index = ((index << 2) & mask) | (basesMap[seq.getBaseAt(i)] & 0x3);
    ((*freq)[index])++;
  }


  uint64_t count = 0;  
  uint64_t unique = 0;
  for (size_t i = 0; i < M; ++i) {
    if ((*freq)[i] > 0) {
      //    std::cout << codeToBases(i,k) << "\t\t" << (*freq)[i] << std::endl;
      count += (*freq)[i];
      unique++;
    }
  }
  std::cout << "Total\t" << count << "\tUnique\t" << unique << std::endl;
  std::cout.flags(f);
  //std::cout << (char*)seq.getSequence() << std::endl;
  

  return freq;
}


