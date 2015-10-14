/**
 * This is a simulator developed ad-hoc to allow future optimization regardless
 * the changes to other components of the library (which are designed ti be part
 * of a library rather then efficient stand alon tools.
 */

#include "common.h"

#include <cstdlib>
#include <ctime>

#include <iostream>
#include <list>

char bases[] = {'A', 'C', 'G', 'T'};
char revBases[128];
Options Options::opts;

void initSimulator() {
  revBases['A'] = revBases['a'] = 0;
  revBases['C'] = revBases['c'] = 1;
  revBases['G'] = revBases['g'] = 2;
  revBases['T'] = revBases['t'] = 3;
  srand(time(NULL));
  initChainMatrix();
}

void clearSimulator() {
  clearChainMatrix();
}

char randomMutation(char c) { 
  return bases[(revBases[c] + ((rand() % 3) + 1) ) & 0x3];
}

void printString(char* s, size_t n) {
  for (int i = 0; i < n; ++i) {
    std::cout << s[i];
  }
}

void simulateReadAt(size_t j, size_t m, char* S, char* r) {
  for (size_t l = 0; l < m; ++l) {
    r[l] = S[j+l];
  }
}

void simpleIIDErrors(std::string& s, double pe) {
  size_t m = s.length();
  for (size_t i = 0; i < m; ++i) {
    double x = (double)rand() / RAND_MAX;
    if (x <= pe) {
      s[i] = randomMutation(s[i]);
    }
  }
}

size_t hammingDistance(const char* s1, const char* s2, size_t m) {
  size_t d = 0;
  for (int i = 0; i < m; ++i) {
    if (s1[i] != s2[i]) {
      ++d;
    }
  }
  return d;
}

size_t hammingDistance(const std::string& s1, const std::string& s2, size_t m) {
  return hammingDistance(s1.c_str(), s2.c_str(), m);
}


int main(int argc, char** argv) {
  std::cout << std::endl;
  
  char* ref = NULL;
  char* read = NULL;
  // Important NOT invert (init requires argument to be parsed)
  parseArguments(argc,argv);
  initSimulator();

  
  size_t N = Options::opts.N;
  size_t m = Options::opts.m;
  size_t M = Options::opts.M;
  size_t Nbar = N - m + 1;

  double pe = 0.1;
  

  std::cout << std::endl;
  std::cout << "\t\t+++++  Starting simulation +++++ \n\n";
  
  std::cout << "* Reference generation... ";
  ref = new char[N];
  generateIIDGenome(N,ref);
  std::cout << "[OK]" << std::endl;

  std::cout << "* Read generation... ";
  std::list<std::string> reads;
  read = new char[m];
  for (size_t h = 0; h < M; ++h) {
    std::string r(ref + (rand() % Nbar),m);
    std::string rr = r;
    simpleIIDErrors(r,pe);
    reads.push_front(r);
  }
  std::cout << "[OK]" << std::endl;
  
  
  std::cout << "* Cleaning... ";

  
  delete[] read;
  delete[] ref;
  //  printChainMatrix();
  clearSimulator();
  std::cout << "[OK]" << std::endl;

  
  std::cout << std::endl;
  return 0;
}
