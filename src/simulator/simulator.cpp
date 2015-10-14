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


int main(int argc, char** argv) {
  std::cout << std::endl;
  
  char* ref = NULL;
  //  char* read = NULL;

  // Important NOT invert (init requires argument to be parsed)
  parseArguments(argc,argv);
  initSimulator();

  
  size_t N = Options::opts.N;
  size_t m = Options::opts.m;
  size_t M = Options::opts.M;
  size_t Nbar = N - m + 1;

  double pe = Options::opts.pe;
  

  std::cout << std::endl;
  std::cout << "\t\t+++++  Starting simulation +++++ \n\n";
  
  std::cout << "* Reference generation... ";
  ref = new char[N];
  generateIIDGenome(N,ref);
  std::string s(ref);
  std::cout << "[OK]" << std::endl;

  std::cout << "* Read generation... ";
  std::list<std::string> reads;
  generateOfflineReads(s, reads);

  
  //  read = new char[m];
  //for (size_t h = 0; h < M; ++h) {
    
    // std::string r(ref + (rand() % Nbar),m);
    // std::string rr = r;
    // simpleIIDErrors(r,pe);
    // reads.push_front(r);
  //  }
  std::cout << "[OK]" << std::endl;
  
  
  std::cout << "* Cleaning... ";

  
  //  delete[] read;
  delete[] ref;
  //  printChainMatrix();
  clearSimulator();
  std::cout << "[OK]" << std::endl;

  
  std::cout << std::endl;
  return 0;
}
