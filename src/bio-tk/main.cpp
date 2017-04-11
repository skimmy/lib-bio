#include "options.hpp"
#include "tasks.hpp"

#include <iostream>

#include <include/prob/generator.hpp>
#include <include/prob/probability.hpp>
#include <generate_task.hpp>

#include <util/io_helper.hpp>

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
  case 5:
    {
      // statistics of reads
      taskSelectedMsg = "Reads statistics";
      string reads = opts.readsFile;
      string w_dir = opts.outputDir;
      string prefix = opts.prefixFile;
      task_read_statistics(reads, w_dir, prefix);
      break;
    }
  case 6:
    {
      using IIDCharSampler = lbio::IIDSampler<lbio::DiscreteProbability<char>>;

      // generation
      // TODO: proper parameters
      lbio_size_t L = 100;
      std::map<char,double> equal_prob_bases =
	{
	  { 'A', 0.25},
	  { 'C', 0.25},
	  { 'G', 0.25},
	  { 'T', 0.25}
	};
      lbio::DiscreteProbability<char,double> _distr(equal_prob_bases.cbegin(),
						    equal_prob_bases.cend());
      IIDCharSampler sampler(_distr);
      std::string seq = generate_iid_bases(L, sampler);
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
  std::stringstream _str_str {"1,Pippo\n2,Plut\n3,John el ton 333"};
  auto _in = lbio::stream_to_map<int,std::string>(_str_str);
  for (auto k : _in) {
    std::cout << k.first << "\t" << k.second << "\n";
  }
}


int main(int argc, char** argv) {
  runTask(argc, argv);
  test();   
  return 0;
}
