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
extern double* power4_lookup;

enum OpMode { Test, Offline, Online, Oracle, ScoreEst };

// options parsing
struct Options {
  size_t N;
  size_t m;
  size_t M;
  double pe;

  size_t empiricalDistributionStep;

  std::string inputReference;
  std::string outputDistribution;
  std::string outputCDF;

  int approxLevel;

  OpMode mode;

  bool pipeline;
  bool online;

  bool verbose;
  bool test;

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
void generateConstantGenome(size_t N, char* S, char b);
void generateOfflineReads(const std::string& s, std::priority_queue<Read>& reads);
Read generateOnlineRead(char* S, size_t j);
void simulateReadAt(size_t j, size_t m, const char* S, char* r);
void simpleIIDErrors(std::string& s, double pe);
char randomMutation(char c);
void complementBases(char* S, size_t n);
char baseComplement(char b);

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
void initUtil();
void clearUtil();

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

// probabilities

// Class representing empirical distribution as histogram array
class EmpiricalDistribution {  
public:
  EmpiricalDistribution(double a, double b, size_t N);

  // given a sample x in [a,b] returns the index where x must be counted
  size_t indexForSample(double x) const;
  // returns the relative frequenxy (i.e., count[i] / total) for the i-th index
  double valueAtIndex(size_t i) const;
  // returns the number of elements at the i-th index
  double countAtIndex(size_t i) const { return f[i]; }

  // returns the number of intervals 'n'
  double getIntervalCount() const { return n; }
  
  // add a new sample x to the distribution
  void addSample(double x);

  // returns the emprical cumulative distribution derived from
  // currently stored samples.
  void getCDF(std::vector<double>& cdf)  const;
    
  
private:
  std::vector<double> f;
  size_t n;
  double xa;
  double xb;
  double step; // defined as (b-a)/n
  double total;
};

struct ScoreSumFreq
{
  double sSum;
  int sFreq;
};

void initProbabilities();
void clearProbabilities();
double randomReadsOverlapProbNoErr(const std::string& s1, const std::string& s2, size_t s);
double overlappingStringsSum(const std::string & s1, const std::string& s2);
double overlappingStringsSumWithErr(const std::string& s1, const std::string& s2);
size_t generateInterReadDistance();

size_t percentileIndex(const std::vector<double>& cdf, double perc);

double score(const std::string& r1, const std::string& r2, size_t s);
double scoreExt(const std::string& r1, const std::string& r2, size_t s, double* num_den);

// testing functions
void testAll();

#endif
