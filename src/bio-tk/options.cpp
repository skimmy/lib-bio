#include "options.hpp"

#include <getopt.h>
#include <cstdlib>

#include <algorithm>
#include <vector>
#include <map>

#include <chrono>
using Clock = std::chrono::system_clock;

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

// TO BE MOVED IN A PROPER TIME/STRING UTIL PLACE

template <typename Num, int Digits>
std::string fixed_len_digit_string(Num n) {
  int d = Digits;
  std::string str(d, '0');
  auto back_str = str.rbegin();
  while(d>0) {
    *back_str = '0' + (n % 10);
    n /= 10;
    ++back_str;
    --d;
  }
  return str;
}

std::string day_from_tm(const std::tm& t) {
  return std::to_string(1900 + t.tm_year) + "-"
    + fixed_len_digit_string<int,2>(1 + t.tm_mon) + "-"
    + fixed_len_digit_string<int,2>(t.tm_mday);
}

std::string time_from_tm(const std::tm& t) {
  return fixed_len_digit_string<int,2>(t.tm_hour) + ":"
    + fixed_len_digit_string<int,2>(t.tm_min)  + ":"
    + fixed_len_digit_string<int,2>(t.tm_sec);
}

string timestamp_string() {
  auto now = Clock::now();
  auto seconds = std::chrono::time_point_cast<std::chrono::seconds>(now);
  auto fraction = now - seconds;
  std::time_t now_time_t = Clock::to_time_t(now);
  std::tm tm = *std::localtime(&now_time_t);
  auto milliseconds =
    std::chrono::duration_cast<std::chrono::milliseconds>(fraction);
  auto msecs = milliseconds.count();
  return day_from_tm(tm) + " " + time_from_tm(tm) + "." + 
    fixed_len_digit_string<decltype(msecs),3>(msecs);
}

size_t
str_position_of(char c, const std::string& text) {
  size_t pos = std::string::npos;
  for(size_t i = 0; i < text.size(); ++i) {
    if (text[i] == c) {
      return i;
    }
  }
  return pos;
}

std::string
remove_punctuation_and_spaces(const std::string& in) {
  std::string banned_chars = ".,;:!? \t-";
  std::string out {};
  auto iBegin = in.cbegin();
  while(iBegin != in.cend()) {
    if (str_position_of(*iBegin, banned_chars) == std::string::npos) {
        out.insert(out.end(), *iBegin);
    }
    ++iBegin;
  }
  return out;
}

// END UTIL FUNCTIONS

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
  prefixFile = remove_punctuation_and_spaces(timestamp_string()) + "_";

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
  if ( argc < 2 ) {
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

  // DEBUG
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
