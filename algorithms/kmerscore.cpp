#include "../algorithms.h"
#include "../sequence.h"

const size_t NoPos = (size_t) -1;

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
    while(i < m) {
      ++firstMap;
      if ( (map[i] != NoPos) && (map[i] != firstMap) ) {
	return false;
      }
      ++i;
    }
  }
  return true;
}
