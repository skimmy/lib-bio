#ifndef COMMON_H
#define COMMON_H

#include <cstdlib>

// general variables
char bases[] = {'A', 'C', 'G', 'T'};
char revBases[128];

// options parsing
struct Options {
  size_t N;
  size_t m;
  size_t M;

  bool verbose;
public:
  static  Options opts;
};
Options Options::opts;
void parseArguments(int argc, char** argv);

// generator
void generateIIDGenome(size_t N, char* S);

#endif
