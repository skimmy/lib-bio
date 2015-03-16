#ifndef GENOME_GENERATOR_H
#define GENOME_GENERATOR_H

#include <memory>


namespace libbio {

  namespace generator {

    std::unique_ptr<char[]> generateIidSequence(size_t N);
    std::unique_ptr<char[]> generateGCRichSequence(size_t N);
    std::unique_ptr<char[]> generateGCPoorSequence(size_t N);
  }
  
}

#endif
