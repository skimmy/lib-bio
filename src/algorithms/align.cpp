#include <vector>
#include <string>

#include "../algorithms.h"

namespace bio {

  std::vector<Position<int>> alignReads(const std::vector<Read>& reads, const std::string& ref) {
    std::vector<Position<int>> aligns;
    size_t n = reads.size();
    for (size_t i = 0; i < n; ++i) {
      Read r = reads[i];
      size_t bestAlign = -1;
      SmithWatermanDP sw(r.getBases().c_str(), r.getBases().size(), ref.c_str(), ref.size());
      sw.computeMatrix();
      bestAlign = sw.getGlobalBest().j;
      aligns.push_back(Position<int>((int)i, bestAlign));
    }
    return aligns;
  }

}
