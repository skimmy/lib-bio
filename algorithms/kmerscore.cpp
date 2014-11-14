#include "../algorithms.h"

const size_t NoPos = (size_t) -1;

KmerScoreType scoreForRead(KmerMap map, size_t k) {
  return 0;
}

bool isKmerUnique(const ReadKmerMap& map) {
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
