#ifndef FASTQ_LAZY_LOADER_H
#define FASTQ_LAZY_LOADER_H

#include <fstream>
#include <list>

#include "FastqRead.hpp"

class FastqLazyLoader {
private:
  std::ifstream input;
public:
  explicit FastqLazyLoader(const string& filePath); 
  ~FastqLazyLoader();
  std::list<FastqRead> getNextReads(const size_t M);
};

#endif
