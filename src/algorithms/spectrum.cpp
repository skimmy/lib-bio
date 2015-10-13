#include "spectrum.hpp"

const size_t NoPos = (size_t) -1;

/**
 */
uint64_t shiftAndPaddWithCharKmer(uint64_t old, char c, size_t k) {
  uint64_t tmp = old << 2;
  tmp &= ( ( 0x1L << (2*k) ) - 1 ) ;
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
  uint64_t n_kmer = seq::NumericKMer::fromChars(kmer, k);
  do {
    index[n_kmer]++;
    n_kmer = shiftAndPaddWithCharKmer(n_kmer, ref.getBaseAt(i), k);
    ++i;
  } while(i <= N);
  delete[] kmer;
  return index;
}

void spectrumAsArray(const Sequence& ref, size_t k, uint64_t* v) {
  // initialize the vector v
  size_t K = 1 << (2 * k);
  memset(v, 0, K * sizeof(uint64_t));
  char* kmer = new char[k];
  // scan the sequence
  size_t N = ref.getSequenceLength();
  size_t i = 0;
  // get the first k-mer
  for (; i < k; ++i) {
    kmer[i] = ref.getBaseAt(i);
  }
  i = k;
  uint64_t n_kmer = seq::NumericKMer::fromChars(kmer, k);
  do {
    v[n_kmer]++;
    n_kmer = shiftAndPaddWithCharKmer(n_kmer, ref.getBaseAt(i), k);
    ++i;
  } while(i <= N);
  delete[] kmer;
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
  uint64_t n_kmer = seq::NumericKMer::fromChars(kmer, k);
  do {
    index[n_kmer].push_back(i);
    n_kmer = shiftAndPaddWithCharKmer(n_kmer, ref.getBaseAt(i), k);
    ++i;
  } while(i <= N);
  delete[] kmer;
  return index;
}

KmersMap extractKmersMapPosition(const Sequence& seq, NumericKmerIndex index, size_t k) {
  size_t m = seq.getSequenceLength();
  size_t barm = m - k + 1;
  KmersMap map = KmersMap(barm);
  char* kmer = new char[k];
  // get the first k-mer
  size_t i = 0;
  for (; i < k; ++i) {
    kmer[i] = seq.getBaseAt(i);
  }
  uint64_t n_kmer = seq::NumericKMer::fromChars(kmer, k);
  for (size_t j = 0; j < barm; ++j) {    
    NumericKmerIndex::const_iterator it = index.find(n_kmer);
    map[j] = (it == index.end()) ? NoPos : (it->second).front();    
    n_kmer = shiftAndPaddWithCharKmer(n_kmer, seq.getBaseAt(j+k), k);
  }
  delete[] kmer;
  return map;
}
