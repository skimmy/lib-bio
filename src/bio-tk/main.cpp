#include "options.hpp"
#include "tasks.hpp"

#include <fstream>

#include <util/io_helper.hpp>

// testing includes
#include <iterator>
#include <algorithms/hamming_distance.hpp>
#include <algorithms/edit_distance.hpp>
#include <../io.h>

#define DEBUG 1

int runTask(OPTIONS& opts) {
  // parse arguments
  

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
      std::ifstream _conf(opts.config_file);
      std::map<std::string, std::string> gen_opts =
	lbio::stream_to_map<std::string, std::string>(_conf, '=');
      for (auto x : gen_opts) {
	std::cerr << x.first << "\t" << x.second << "\n";
      }
      task_generate(gen_opts);
      break;
    }
  default:
    std::cout << "Unrecognized operation" << std::endl;
    return 1;
  }
  std::cout << "[TASK]    " << taskSelectedMsg << std::endl;
  return 0;
} 

void test(OPTIONS& opts) {
  // Test disclaimer
  std::cout << "\n********** WARNING ********** \n" <<  
    "  This is a Test Release \n" << 
    "***************************** \n" <<  std::endl;
  
  std::ifstream _fastq1("/home/skimmy/data/EDAF/iid_test_1.fastq");
  std::ifstream _fastq2("/home/skimmy/data/EDAF/iid_test_2.fastq");
  std::istream_iterator<FastqRead> _it1(_fastq1);
  std::istream_iterator<FastqRead> _it2(_fastq2);
  std::istream_iterator<FastqRead> _it_end;
  lbio::edit_distance_wf<std::string> edit(500,500);
  while ( (_it1 != _it_end) and (_it2 != _it_end) ) {
    FastqRead r1 {*_it1};
    FastqRead r2 {*_it2};
    std::string s1 = r1.getBases();
    std::string s2 = r2.getBases();
    lbio_size_t d = lbio::hamming_distance(s1.cbegin(), s1.cend(), s2.cbegin());
    lbio_size_t e = edit.compute(s1, s1.size(), s2, s2.size());
    std::cout << s1 << "\n" << s2 << "\n" << d << "\t" << e << "\n";
    ++_it1;
    ++_it2;
  }
}


int main(int argc, char** argv) {
  OPTIONS opts;
  opts.parseInputArgs(argc, argv);
#ifndef DEBUG
  opts.printOptions(std::cout);
#endif
  runTask(opts);
  test(opts);   
  return 0;
}
