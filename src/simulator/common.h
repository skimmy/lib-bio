#ifndef COMMON_H
#define COMMON_H

#include <cstdlib>
#include <string>
#include <queue>

// general variables
extern char bases[];
extern char revBases[128];

// options parsing
struct Options {
  size_t N;
  size_t m;
  size_t M;

  double pe;

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
  
void generateIIDGenome(size_t N, char* S);
void generateOfflineReads(const std::string& s, std::priority_queue<Read>& reads);
void simulateReadAt(size_t j, size_t m, const char* S, char* r);
void simpleIIDErrors(std::string& s, double pe);
char randomMutation(char c);

// chain
void initFalsePositiveMatrix();
void clearFalsePositiveMatrix();
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

void printString(char* s, size_t n);
 
void printDoubleMatrix(double** M, size_t n, size_t m);
double** initDoubleMatrix(size_t n, size_t m);
void clearDoubleMatrix(double** matrix, size_t n, size_t m);
double elementsSumDoubleMatrix(double** matrix, size_t n, size_t m);

// probabilities
void initProbabilities();
void clearProbabilities();

#endif
