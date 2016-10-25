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



  
}

#include "../core/DNAAlphabet2Bits.hpp"

int main(int argc, char** argv) {
  //runTask(argc, argv);
  test();   
  return 0;
}
