#include <memory>
#include <random>

namespace libbio {
  namespace generator {

    std::discrete_distribution<int> iidDist {1,1,1,1};
    std::discrete_distribution<int> GCRichDist {1,2,2,1};
    std::discrete_distribution<int> GCPoorDist {2,1,1,2};
    std::default_random_engine generator;

    std::unique_ptr<char[]>  generateSequenceWithDistribution
    (size_t N, std::discrete_distribution<int> & dist) {
      std::unique_ptr<char[]> pSeq(new char[N]);    
      char bases[] = { 'A', 'C', 'G', 'T' };
      for (size_t i = 0; i < N; ++i) {
	pSeq[i] = bases[dist(generator)];
      }
      return pSeq;
    }

    std::unique_ptr<char[]> generateIidSequence(size_t N) {    
      return generateSequenceWithDistribution(N, iidDist);
    }

    std::unique_ptr<char[]> generateGCRichSequence(size_t N) {    
      return generateSequenceWithDistribution(N, GCRichDist);
    }

    std::unique_ptr<char[]> generateGCPoorSequence(size_t N) {    
      return generateSequenceWithDistribution(N, GCPoorDist);
    }

  }
}
