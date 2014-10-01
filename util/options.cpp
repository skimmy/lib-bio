#include "options.hpp"

#include <getopt.h>
#include <stdlib.h>

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

const char* shortOptions = "hntg:G:f:r:R:F:o:p:c:A:k:";
const struct option longOptions[] = 
  {
    { "help", 0, NULL, 'h' },
    { "no-align", 0, NULL, 'n' },
    { "translate", 0, NULL, 't' },
    { "genome", 1, NULL, 'g' },
    { "output-genome", 1, NULL, 'G' },
    { "genome-format", 1, NULL, 'f' },
    { "reads", 1, NULL, 'r' },
    { "output-reads", 1, NULL, 'R' },
    { "reads-format", 1, NULL, 'F' },
    { "output-align", 1, NULL, 'o' },
    { "padding", 1, NULL, 'p' },
    { "genome-copies", 1, NULL, 'c' },
    { "algorithm-type", 1, NULL, 'A' },
    { "kmer-size", 1, NULL, 'k' }
  };

options::options() {
  genomeFormat = GENOME_CUSTOM;
  readsFormat = READS_CUSTOM;
  genomeFile = "genome.dat";
  readsFile = "reads.dat";
  
  genomeOutputFile = "genome.dat";
  readsOutputFile = "reads.dat";
  alignOutputFile = "align.out";

  padding = 0;
  genomeCopies = 1;
  
  translate = false;
  align = true;
  alignAlgorithm = CPU_DP;
  kmerSize = 15;
}

void options::printUsage(ostream& os, const char* name, int exitCode) {
  os << "Usage:\n\t" << name << " [OPTIONS] [FILE(S)]\n";
  exit(exitCode);
}

void options::parseInputArgs(int argc, char** argv) {
  if ( argc < 3 ) {
    this->printUsage(cerr, argv[0], 1);
  }
  int nextOption = -2;
  do {
    nextOption = getopt_long( argc, argv, shortOptions, longOptions, NULL );
    switch(nextOption) {
    case 'h':
      printUsage(cout, argv[0], 0);
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
    default:
      return;
      //      this->printUsage(cerr, argv[0], 1);
    }
  } while( nextOption != -1);
}

void options::printOptions(ostream& os) {
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
}
