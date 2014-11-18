#include "../algorithms.h"
#include "../sequence.h"

const size_t NoPos = (size_t) -1;

// Get a reference and a read and compute the vector
//   A = a1 a2 .. 
// where ai = 1 iff i-th kmer maps (somewhere)

// Get index for a reference and compute the mapping
//   P = p1 p2 ...
// where pi is s.t. Read[i,...,i+k - ] = Ref[pi, ...] or pi = -1

bool isKmerUnique(const KmersMap& map) {
  if (map.empty()) {
    return true;
  }
  size_t n = map.size();
  size_t prev = map[0];
  for (size_t i = 1; i < n; ++i) {
    if(map[i] != NoPos && map[i] != (prev + 1)) {
      return false;
    }
  }
  return true;
}
