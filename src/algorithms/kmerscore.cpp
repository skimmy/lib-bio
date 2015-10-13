#include "../algorithms.h"
#include "../sequence.h"

#include <algorithm>

const size_t NoPos = (size_t) -1;

void inline add_k(std::vector< uint64_t >& v, size_t i, size_t k, uint64_t o = 1) {
  size_t n = std::min(v.size(), i + k);
  while(i < n) {
    v[i] += o;
    ++i;
  }
}

std::vector< uint64_t > kmerScoreVector(const KmersMap& map, size_t k) {
  size_t barm = map.size();
  size_t m = barm + k - 1;
  std::vector< uint64_t > v(m, 0);
  for(size_t i = 0; i < barm; ++i) {
    if (map[i] != NoPos) {
      add_k(v, i, k);
    }
  }
  return v;
}

KmerScoreType scoreForVector(const std::vector< uint64_t >& v, size_t k) {
  KmerScoreType score = 0;
  size_t barm = v.size() - k + 1;
  std::vector< uint64_t >::const_iterator it = v.begin();
  while(it != v.end()) {
    score += (KmerScoreType)(*it);
    ++it;
  }
  return score / ((double)(k*barm));
}

size_t kmerErrorCount(const KmersMap& map, size_t k) {
  std::vector< uint64_t > v = kmerScoreVector(map, k);
  size_t errs = 0;
  std::vector< uint64_t >::const_iterator it = v.begin();
  while(it != v.end()) {    
    errs += (size_t)((*it) == 0);
    ++it;
  }
  return errs;
}

// Get a reference and a read and compute the vector
//   A = a1 a2 .. 
// where ai = 1 iff i-th kmer maps (somewhere)

// Get index for a reference and compute the mapping
//   P = p1 p2 ...
// where pi is s.t. Read[i,...,i+k - ] = Ref[pi, ...] or pi = -1

bool isKmerUniquelyMapped(const KmersMap& map) {
  if (map.empty()) {
    return true;
  }
  size_t firstMap = 0;
  size_t m = map.size();
  size_t i = 0;
  // search the first position with map not equal to -1
  while(map[i] == NoPos && i < m) {
    ++i;
  }
  if (i < m) {
    firstMap = map[i];
    ++i;
    // scan till the and
    while(i < m) {
      ++firstMap;
      // if map[i] is different from -1 it must match the current counter
      if ( (map[i] != NoPos) && (map[i] != firstMap) ) {
	return false;
      }
      ++i;
    }
  }
  return true;
}
