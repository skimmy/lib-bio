#include "options.hpp"

#include <getopt.h>
#include <cstdlib>

#include <algorithm>
#include <vector>
#include <map>

const std::vector<std::string> tVec = 
  {
    std::string("nop"),
    std::string("align"),
    std::string("kspectrum"),
    std::string("kmapping"),
    std::string("kscore"),
    std::string("readstats")
  };


int parseTask(const string& taskName) {
  auto pos = std::find(tVec.begin(), tVec.end(), taskName);
  return (pos != tVec.end()) ?
    std::distance(tVec.begin(),pos) : -1;
}

GenomeFormat parseGenomeFormat(const char* optval) {
  int v = atoi(optval);
  switch(v) {
  case 0:
    return GENOME_CUSTOM;
  case 1:
    return FAST;
  default:
    return GENOME_CUSTOM;
  }
}

ReadsFormat parseReadsFormat(const char* optval) {
  int v = atoi(optval);
  switch(v) {
  case 0:
    return READS_CUSTOM;
  case 1:
    return FASTQ;
  case 2:
    return CSFASTA;
  default:
    return READS_CUSTOM;
  }
}

AlignAlgorithm parseAlgorithmType(const char* optval) {
  int v = atoi(optval);
  switch(v) {
  case 0:
    return CPU_DP;
  case 1:
    return GPU_DP;
  default:
    return CPU_DP;
  }
}

string timestamp_string() {
  return "";
}

const char* shortOptions = "hvntg:G:f:r:R:F:o:d:X:p:c:A:k:T:";
const struct option longOptions[] = 
  {
    { "help", 0, NULL, 'h' },
    { "verbose", 0, NULL, 'v' },
    { "no-align", 0, NULL, 'n' },
    { "translate", 0, NULL, 't' },
    { "genome", 1, NULL, 'g' },
    { "output-genome", 1, NULL, 'G' },
    { "genome-format", 1, NULL, 'f' },
    { "reads", 1, NULL, 'r' },
    { "output-reads", 1, NULL, 'R' },
    { "reads-format", 1, NULL, 'F' },
    { "output-align", 1, NULL, 'o' },
    { "output-dir", 1, NULL, 'd' },
    { "file-name-prefix", 1, NULL, 'X' },
    { "padding", 1, NULL, 'p' },
    { "genome-copies", 1, NULL, 'c' },
    { "algorithm-type", 1, NULL, 'A' },
    { "kmer-size", 1, NULL, 'k' },
    { "threads", 1, NULL, 'T' }
  };

options::options() {

  task = 0;
  verbose = false;

  genomeFormat = GENOME_CUSTOM;
  readsFormat = READS_CUSTOM;
  genomeFile = "genome.dat";
  readsFile = "reads.dat";
  
  genomeOutputFile = "genome.dat";
  readsOutputFile = "reads.dat";
  alignOutputFile = "";

  outputDir = "/tmp/";
  prefixFile = timestamp_string() + "_";

  padding = 0;
  genomeCopies = 1;

  
  translate = false;
  align = true;
  alignAlgorithm = CPU_DP;
  kmerSize = 15;

  threadsNumber = 1;

  
}

void options::printUsage(ostream& os, const char* name, int exitCode) {
  os << "Usage:\n\t" << name << " [OPTIONS] [FILE(S)]\n";
  exit(exitCode);
}

void options::parseInputArgs(int argc, char** argv) {
  if ( argc < 4 ) {
    this->printUsage(cerr, argv[0], 1);
  }
  task = parseTask(argv[1]);
  int nextOption = -2;
  do {
    nextOption = getopt_long( argc, argv, shortOptions, longOptions, NULL );
    switch(nextOption) {
    case 'h':
      printUsage(cout, argv[0], 0);
      break;
    case 'v':
      this->verbose = true;
      break;
    case 'n':
      this->align = false;
      break;
    case 't':
      this->translate = true;
      break;
    case 'g':
      this->genomeFile = optarg;
      break;
    case 'G':
      this->genomeOutputFile = optarg;
      break;
    case 'f':
      this->genomeFormat = parseGenomeFormat(optarg);
      break;
    case 'r':
      this->readsFile = optarg;
      break;
    case 'R':
      this->readsOutputFile = optarg;
      break;
    case 'F':
      this->readsFormat = parseReadsFormat(optarg);
    case 'o':
      this->alignOutputFile = optarg;
      break;
    case 'd':
      this->outputDir = optarg;
      break;
    case 'X':
      this->prefixFile = optarg;
      break;
    case 'p':
      this->padding = atoi(optarg);
    case 'c':
      this->genomeCopies = atoi(optarg);
    case 'A':
      this->alignAlgorithm = parseAlgorithmType(optarg);
      break;
    case 'k':
      this->kmerSize = atoi(optarg);
      break;
    case 'T':
      this->threadsNumber = atoi(optarg);
    default:
      return;
      //      this->printUsage(cerr, argv[0], 1);
    }
  } while( nextOption != -1);
}

void options::printOptions(ostream& os) {
  os << "-- TASK \n";
  os << tVec[this->task] << '\n';
  os << "-- INPUT --\n";
  os << "Genome: " << this->genomeFile << '\n';
  os << "Reads:  " << this->readsFile << '\n';
  os << "-- OUTPUT --\n";
  os << "Genome: " << this->genomeOutputFile << '\n';
  os << "Reads:  " << this->readsOutputFile << '\n';
  os << "Align:  " << this->alignOutputFile << '\n';
  os << "-- PRE PROCESSING\n";
  os << "Padding:       " << this->padding << '\n';
  os << "Genome copies: " << this->genomeCopies << '\n';
  os << "-- ALGORITHM --\n";
  os << "Translate:  " << this->translate << '\n';
  os << "Alignment:  " << this->align << '\n';
  os << "Algorithm:  " << this->alignAlgorithm << '\n';
  os << "K-Mer size: " << this->kmerSize << '\n';
  os << "-- RUNTIME --" << '\n';
  os << "Threads:    " << this->threadsNumber;
}
