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
  initProbabilities();
  initChainMatrix();
  initFalsePositiveMatrix();
}

void clearSimulator() {
  clearFalsePositiveMatrix();
  clearChainMatrix();
  clearProbabilities();
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
  double p_fail = 1.0;
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
      size_t est_s = bestHammingOverlap(r2.r, r1.r);
      size_t obs_dh = prefixSuffixHammingDistance(r1.r, r2.r, est_s);
      //      p_fail += 1.0 - randomReadsOverlapProbNoErr(est_s,obs_dh);
      p_fail += 1.0 - randomReadsOverlapProbNoErr(s,dh);
    } else {
      addNonOverlapRecord(r2.j - r1.j - m);
    }
    r1 = r2;
  }
  std::cout << "[OK]" << std::endl;

  std::cout << "P[Fail]    = " << p_fail << std::endl;
  std::cout << "P[Success] = " << 1.0 - p_fail << std::endl;
  
  std::cout << "* Cleaning... ";
  delete[] ref;
  std::cout << "[OK]" << std::endl;

  
  std::cout << std::endl;
  return 0;
}
