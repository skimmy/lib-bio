#include "options.hpp"
#include "tasks.hpp"
#include "edaf_task.hpp"

#include <fstream>

#include <util/io_helper.hpp>

// testing includes
//#include <util/str_util.hpp>
#include <generate_task.hpp>

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
  case 6: // generate
    {
      // data generation
      std::ifstream _conf(opts.config_file);
      std::map<std::string, std::string> gen_opts =
	lbio::stream_to_map<std::string, std::string>(_conf, '=');
      for (auto x : gen_opts) {
	std::cerr << x.first << "\t" << x.second << "\n";
      }
      task_generate(gen_opts);
      break;
    }
  case 7: // edaf (Edit Distance VS Alignment Free)
    {
      std::cout << "EDAF\n";
      std::ifstream _conf(opts.config_file);
      std::map<std::string, std::string> gen_opts =
	lbio::stream_to_map<std::string, std::string>(_conf, '=');
      edaf_task(gen_opts);
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
  std::vector<double> p_ops {0.06, 0.03, 0.03}; // Sub ~6% Del=Ins ~3% -> Mat ~88%
  std::string read = "ACCGTTGCGTAAAGAT";
  std::string noisy_read = generate_iid_edit_transformation(read, p_ops);
  std::cout << std::endl << read << std::endl;
  std::cout << noisy_read << std::endl;
  
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
