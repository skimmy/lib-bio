#ifndef COMMON_H
#define COMMON_H

#include <cstdlib>
#include <string>
#include <queue>
#include <vector>
#include <fstream>

// general variables
extern char bases[];
extern char revBases[128];

// options parsing
struct Options {
  size_t N;
  size_t m;
  size_t M;

  double pe;

  std::string inputReference;  

  bool pipeline;
  bool online;

  bool verbose;  

  static  Options opts;
};

void parseArguments(int argc, char** argv);

// generator
struct Read {
  std::string r;
  size_t j;
  Read(const std::string& s, size_t p) : j(p), r(s) {}
  // !!! WARNING '>' is used to make the priority queue work in ascending order
  bool operator < (const Read& other) const { return this->j > other.j; }
};

// used to represent currenct 'under sequence' segment of a large genome (used
// when simulator is operating in 'online' mode)
struct GenomeSegment {
  char* genome;

  size_t length;
  size_t current_start;
  size_t total_length;
  size_t read_length;
  
  
  GenomeSegment(size_t N, size_t m, size_t l);
  ~GenomeSegment();
};
  
void generateIIDGenome(size_t N, char* S);
void generateOfflineReads(const std::string& s, std::priority_queue<Read>& reads);
Read generateOnlineRead(char* S, size_t j);
void simulateReadAt(size_t j, size_t m, const char* S, char* r);
void simpleIIDErrors(std::string& s, double pe);
char randomMutation(char c);

const size_t MAX_GENOME_SEGMENT_LENGTH = 1 << 20;
void generateFirstGenomeSegment(GenomeSegment& g);
void generateNewGenomeSegment(GenomeSegment& g, size_t keep_pref);

// chain
void printFalsePositiveMatrix();

void initChainMatrix();
void clearChainMatrix();
void printChainMatrix();
void printNonOverlapDistribution();
void evaluateChainRelation(const Read& r1, const Read& r2, size_t s);
void addNonOverlapRecord(size_t d);

void computeExpectedFalsePositiveMatrix(double** m);

// util
size_t hammingDistance(const char* s1, const char* s2, size_t m);
size_t hammingDistance(const std::string& s1, const std::string& s2, size_t m);
size_t prefixSuffixHammingDistance(const std::string& s1, const std::string& s2, size_t k);
size_t bestHammingOverlap(const std::string& s1, const std::string& s2);

void printString(char* s, size_t n);
void printDoubleMatrix(double** M, size_t n, size_t m);
double** initDoubleMatrix(size_t n, size_t m);
void clearDoubleMatrix(double** matrix, size_t n, size_t m);
double elementsSumDoubleMatrix(double** matrix, size_t n, size_t m);

Read randomRead(size_t m);

std::ifstream openFastaFile(const std::string & path);

// probabilities

class EmpiricalDistribution {
public:
  EmpiricalDistribution(double a, double b, size_t N);
  
private:
  std::vector<double> f;
  size_t n;
  double xa;
  double xb;
};

void initProbabilities();
void clearProbabilities();
double randomReadsOverlapProbNoErr(const std::string& s1, const std::string& s2, size_t s);
double overlappingStringsSum(const std::string & s1, const std::string& s2);
size_t generateInterReadDistance();

#endif
