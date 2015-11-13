/**
 * This is a simulator developed ad-hoc to allow future optimization regardless
 * the changes to other components of the library (which are designed ti be part
 * of a library rather then efficient stand alon tools.
 */

#include "common.h"

#include <cstdlib>
#include <ctime>

#include <iostream>


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
  initFalsePositiveMatrix();
}

void clearSimulator() {
  clearChainMatrix();
  clearFalsePositiveMatrix();
}


int main(int argc, char** argv) {
  std::cout << std::endl;
  
  char* ref = NULL;

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
  std::priority_queue<Read> reads;
  generateOfflineReads(s, reads);
  std::cout << "[OK]" << std::endl;

  std::cout << "* Processing reads... ";

  Read r1 = reads.top();
  reads.pop();
  while(!reads.empty()) {
    Read r2 = reads.top();
    reads.pop();
    size_t s = m - (r2.j - r1.j);
    evaluateChainRelation(r1, r2, s);
    int dh = -1;
    if (s <= m) {
      dh = prefixSuffixHammingDistance(r2.r, r1.r, s);
    } else {
      addNonOverlapRecord(r2.j - r1.j - m);
    }
    r1 = r2;
  }

  std::cout << "[OK]" << std::endl;

  
  std::cout << "* Cleaning... ";
  delete[] ref;
  std::cout << "\n\n";
  printChainMatrix();
  std::cout << "\n\n";
  //printNonOverlapDistribution();
  printFalsePositiveMatrix();
  std::cout << "\n\n";
  double** exp = initDoubleMatrix(m, m+1);
  computeExpectedFalsePositiveMatrix(exp);
  printDoubleMatrix(exp, m, m+1);
  double sum = elementsSumDoubleMatrix(exp, m, m+1);
  std::cout << "P[FP] = " << sum << '\n';
  clearDoubleMatrix(exp, m,m+1);  
  clearSimulator();
  std::cout << "[OK]" << std::endl;

  
  std::cout << std::endl;
  return 0;
}
