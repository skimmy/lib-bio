#include "spectrum.hpp"

/**
 */
uint64_t shiftAndPaddWithCharKmer(uint64_t old, char c, size_t k) {
  uint64_t tmp = old << 2;
  tmp &= ( ( 0x1 << (2*k) ) - 1 ) ;
  tmp |= (DNAAlphabet2Bits::charToInt(c) & 0x3) ;
  return tmp;  
}

std::unordered_map< uint64_t, uint64_t > spectrumAsIntMap(const Sequence& ref, size_t k) {
  // construct the map
  unordered_map< uint64_t, uint64_t > index;
  char* kmer = new char[k];
  // scan the sequence
  size_t N = ref.getSequenceLength();
  size_t i = 0;
  // get the first k-mer
  for (; i < k; ++i) {
    kmer[i] = ref.getBaseAt(i);
  }
  i = k;
  uint64_t n_kmer = NumericKMer::fromChars(kmer, k);
  do {
    index[n_kmer]++;
    n_kmer = shiftAndPaddWithCharKmer(n_kmer, ref.getBaseAt(i), k);
    ++i;
  } while(i <= N);
  delete[] kmer;
  return index;
}


std::unordered_map< uint64_t, std::list< size_t > > kmersMapping(const Sequence& ref, size_t k) {
  unordered_map< uint64_t, std::list< size_t > > index;
  char* kmer = new char[k];
  // scan the sequence
  size_t N = ref.getSequenceLength();
  size_t i = 0;
  // get the first k-mer
  for (; i < k; ++i) {
    kmer[i] = ref.getBaseAt(i);
  }
  i = k;
  uint64_t n_kmer = NumericKMer::fromChars(kmer, k);
  do {
    index[n_kmer].push_back(i);
    n_kmer = shiftAndPaddWithCharKmer(n_kmer, ref.getBaseAt(i), k);
    ++i;
  } while(i <= N);
  delete[] kmer;
  return index;
}
