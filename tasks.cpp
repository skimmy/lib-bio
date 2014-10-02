#include <fstream>

#include "tasks.hpp"


std::vector<Position<int>> alignFastqReadsSimpleSW(const string& readsPath, const string& referencePath, 
						   uint nThreads, size_t nReads) {
  std::vector<Position<int>> aligns;
  // some constants
  uint T = (nThreads > 0) ? nThreads : 1 ;  

  // open files
  std::ifstream readsIn(readsPath, std::ifstream::in);
  std::ifstream refIn(referencePath, std::ifstream::in);

  // read input reference...

  // ...and reads

  // each thread can start as soon as right amount of reads are available
  return aligns;
}
