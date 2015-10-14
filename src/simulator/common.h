#ifndef COMMON_H
#define COMMON_H

#include <cstdlib>
#include <string>
#include <list>

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
void generateIIDGenome(size_t N, char* S);
void generateOfflineReads(const std::string& s, std::list<std::string>& reads);
void simulateReadAt(size_t j, size_t m, const char* S, char* r);
void simpleIIDErrors(std::string& s, double pe);
char randomMutation(char c);

// chain
void initChainMatrix();
void clearChainMatrix();
void printChainMatrix();

// utili
size_t hammingDistance(const char* s1, const char* s2, size_t m);
size_t hammingDistance(const std::string& s1, const std::string& s2, size_t m);

void printString(char* s, size_t n);

#endif
