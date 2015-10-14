#ifndef COMMON_H
#define COMMON_H

#include <cstdlib>

// general variables
extern char bases[];
extern char revBases[128];

// options parsing
struct Options {
  size_t N;
  size_t m;
  size_t M;

  bool verbose;
public:
  static  Options opts;
};

void parseArguments(int argc, char** argv);

// generator
void generateIIDGenome(size_t N, char* S);

// chain
void initChainMatrix();
void clearChainMatrix();
void printChainMatrix();

#endif
