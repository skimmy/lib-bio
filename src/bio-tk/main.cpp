/*#include "sequence.h"
#include "io.h"
#include "util.h"
#include "alignment.h"
#include "quality.h"
#include "adt.h"
#include "tasks.hpp"
#include "algorithms.h"
//#include "generator.h"
//#include "filtering.h"
*/

#include "options.hpp"
#include "tasks.hpp"

#include "../core.h"
#include "../io.h"
#include "../algorithms.h"


#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unordered_map>
#include <queue>

#include <cstdlib>
#include <cmath>

#define DEBUG 1

int runTask(int argc, char** argv) {
  // parse arguments
  OPTIONS opts;
  opts.parseInputArgs(argc, argv);
#ifndef DEBUG
  opts.printOptions(std::cout);
#endif
  int task = opts.task;
  std::string taskSelectedMsg = "None";
  switch(task) {
  case 0:
    {      
      taskSelectedMsg = "No operation";
      break;
    }
  case 1: // Smith Waterman alignment
    {
      taskSelectedMsg = "Smith-Waterman alignment";
      string ref = opts.genomeFile; 
      string reads = opts.readsFile; 
      alignFastqReadsSimpleSW(reads, ref, std::cout, 2, 8);      
      break;
    }
  case 2: // k-spectrum calculation
    {
      taskSelectedMsg = "k-spectrum";
      size_t k = opts.kmerSize;
      string ref = opts.genomeFile;
      taskComputeKSpectrum(k, ref);
      break;
    }
  case 3:
    {
      
      // k-mer mapping
      taskSelectedMsg = "k-mapping";
      size_t k = opts.kmerSize;
      string ref = opts.genomeFile;
      string reads = opts.readsFile;
      string out = opts.alignOutputFile;
      taskMapReadsKmers(ref, reads, k, out);
      break;
    }
  case 4:
    {
      // k-mer score for reads
      taskSelectedMsg = "k-score";
      size_t k = opts.kmerSize;
      string ref = opts.genomeFile;
      string reads = opts.readsFile;
      string out = opts.alignOutputFile;
      size_t nThreads = opts.threadsNumber;
      taskKmerScoreReads(ref, reads, k, out, nThreads);
      break;
    }
  default:
    std::cout << "Unrecognized operation" << std::endl;
    return 1;
  }
  std::cout << "[TASK]    " << taskSelectedMsg << std::endl;
  return 0;
} 

void test() {
  // Test disclaimer
  std::cout << "\n********** WARNING ********** \n" <<  
    "  This is a Test Release \n" << 
    "***************************** \n" <<  std::endl;

  size_t k = 13;
  FastFormat fast("/home/skimmy/biodft/data/seq_iid.fasta");
  Reference s = fast.toReference();
  std::list<KMer> kmers = s.getKMerList(k);


  std::priority_queue< std::pair< double, size_t > > dftVsHamming;
  std::vector< double > dftMods(kmers.size());
  
  size_t i = 0, j = 0;
  double epsilon = 0.5;
  std::vector< uint64_t > hammDist(k);

  // pre-compute all the modules for dft components
  for (std::list<KMer>::iterator it = kmers.begin(); it != kmers.end(); it++) {
    double mod = std::abs(dft::dftComponent(dft::basesToComplexVector(*it), 1, k ));    
    dftMods[i++] = mod;
  }


  // all VS all N^2 scan
  i = j = 0;
  //  size_t MAX_Q = 1 << 22;
  for (std::list<KMer>::iterator it = kmers.begin(); std::next(it) != kmers.end(); it++) {    
    j = i+1;
    for (std::list<KMer>::iterator it2 = std::next(it); it2 != kmers.end(); it2++) {
      
      double d_dft = std::abs(dftMods[i] - dftMods[j]);
      size_t d_h = bio::hammingDistance(*it, *it2);
      if (d_dft < epsilon) {
	hammDist[d_h]++;
      }
      j++;
    }
    i++;
  }

  for (i = 0; i < k; ++i) {
    std::cout << i << " " << hammDist[i] << std::endl;
  }

  
}

#include "../core/DNAAlphabet2Bits.hpp"

int main(int argc, char** argv) {
  //runTask(argc, argv);
  test();   
  return 0;
}
