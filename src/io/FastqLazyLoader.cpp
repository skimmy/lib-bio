#include "../io.h"

FastqLazyLoader::FastqLazyLoader(const string& filePath) 
  : input(filePath) {
  
}

FastqLazyLoader::~FastqLazyLoader() {
  input.close();
}

std::list<FastqRead> FastqLazyLoader::getNextReads(const size_t M) {
  std::list<FastqRead> theList;
  size_t i = 0;
  while (i < M && !(this->input.eof())) {
    FastqRead newRead;
    this->input >> newRead;
    theList.push_back(newRead);
    ++i;
  }
  return theList;
}
