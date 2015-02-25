#include "../io.h"

FastqLazyLoader::FastqLazyLoader(const string& filePath) 
  : input(filePath) {
  
}

FastqLazyLoader::~FastqLazyLoader() {
  input.close();
}

std::list<FastqRead> FastqLazyLoader::getNextReads(const size_t M) {
  std::list<FastqRead> theList;

  return theList;
}
