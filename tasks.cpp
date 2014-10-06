#include <math.h>

#include <vector>
#include <fstream>
#include <thread>

#include "tasks.hpp"
#include "sequence.h"
#include "io.h"
#include "alignment.h"

// This utilitu function returns the bytes remaining until the end of the file
// is reached. It should be probably made publicly availale in a separate util
// section of the library (e.g. util/io.h)
size_t bytesRemainingToTheEnd(ifstream& ifs) {
  std::streampos tmp = ifs.tellg();
  ifs.seekg(0, ifs.end);
  size_t totalLength = ifs.tellg();
  ifs.seekg(tmp);
  return (totalLength - ifs.tellg());
}

void alignSmithWaterman(std::vector<Read>* reads, const Reference* ref) {
  int nReads = reads->size();
  string refBases((char*)ref->getSequence());
  std::cout << "Aligning " << nReads << " reads" << std::endl;
  for (int i = 0; i < nReads; ++i) {    
    SmithWatermanDP sw((*reads)[i].getBases(), refBases);
    sw.computeMatrix();
    MatrixPoint2D maxP = sw.getGlobalBest();
    std::cout << i << " - "  << maxP.i << ", " << maxP.j << std::endl;
  }
}


std::vector<Position<int>> alignFastqReadsSimpleSW(const string& readsPath, const string& referencePath, 
						   uint64_t nThreads, size_t nReads) {

  std::cout << "-- Smith Waterman alignment --" << std::endl;
  std::cout.flush();
  
  std::vector<Position<int>> aligns;
  // some constants
  uint64_t T = (nThreads > 0) ? nThreads : 1;
  uint64_t M = (nReads > 0) ? nReads : 0;

  // open files
  std::ifstream readsIn(readsPath, std::ifstream::in);
  //  std::ifstream refIn(referencePath, std::ifstream::in);  

  // list of all threads (use later for joining)  
  std::vector<std::thread> threads;

  // read input reference...
  std::cout << "    Loading reference..." << std::endl;
  std::cout.flush();
  FastFormat fast;
  fast.loadFromFile(referencePath);
  Reference ref = (Reference)fast;

  

  // ...and reads
  std::cout << "    Loading reads..." << std::endl;
  std::cout.flush();
  if (nReads >= 0) {
    // CASE 1: A specific total number of reads has been specified or
    uint64_t readsPerThread = (size_t)ceil( ((double)(M)) / ((double)(T)) );
    uint64_t actualReads = 0;    
    // repeat for each thread
    for (uint64_t t = 0; t < T; ++t) {
      std::vector<Read>* reads = new std::vector<Read>();
      actualReads = 0;
      while( (!readsIn.eof()) &&  (actualReads < readsPerThread) ) {	
	FastqRead r;
	readsIn >> r;
	reads->push_back(r);
      	actualReads++;
      }
      threads.push_back(std::thread(alignSmithWaterman,reads,&ref));
    }


  } else {
    // CASE 2: we need to estimate the amount of reads to be assigned to each thread
    double fraction = 0;
    double D = 1.0 / ((double)T);
    // repeat for each thread
    for (uint64_t t = 0; t < T; ++t) {
      while( (!readsIn.eof()) && (fraction < D)) {
      	fraction += 0.1;
      }
      //      threads.push_back(std::thread(alignSmithWaterman,reads,&ref));
    }
    
  }

  // synchronize all threads
  for (auto& thr : threads) {
    thr.join();
  }
  return aligns;
}
