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

void offlineSimulation() {

  char* ref = NULL;
  
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

  // Temporary variables to count the number of holes, in the future a more
  // sophisticated way (e.g., finite state machine) should be used.
  size_t holes = 0;
  bool onHole = false;
  
  double p_fail = 0;
  Read r1 = reads.top();
  reads.pop();
  while(!reads.empty()) {
    Read r2 = reads.top();
    reads.pop();
    size_t s = m - (r2.j - r1.j);
    evaluateChainRelation(r1, r2, s);
    int dh = -1;

    if (s <= m) {
      onHole = false;
      double p_ab = randomReadsOverlapProbNoErr(r1.r,r2.r,s);
      if (p_ab < 0.1) {
	std::cout << "\n" << p_ab << " (" << s << ")\t\t" << r1.r << '\t' << r2.r << '\n';
      }
      p_fail += 1.0 - p_ab;

    } else {
      if (onHole == false) {
	onHole = true;
	holes++;
      }
      addNonOverlapRecord(r2.j - r1.j - m);
      double x = (double)Options::opts.N - 2.0 * (double)Options::opts.m + 1.0
	+ overlappingStringsSum(r1.r, r2.r);
      p_fail += 1.0 - ( 1.0 / x);
    }
    r1 = r2;
  }
  std::cout << "[OK]" << std::endl;

  std::cout << "P[Fail]    = " << p_fail << std::endl;
  std::cout << "P[Success] = " << 1.0 - p_fail << std::endl;

  std::cout << "#[Holes]   = " << holes << std::endl;
  
  std::cout << "* Cleaning... ";
  delete[] ref;
  std::cout << "[OK]" << std::endl;

}

void onlineSimulation() {  
  const size_t MAX_GENOME_LENGTH = 2 << 20;

  size_t N = Options::opts.N;
  size_t m = Options::opts.m;
  size_t M = Options::opts.M;
  
  //  GenomeSegment g(N, m, MAX_GENOME_LENGTH);
  GenomeSegment g(N, m, 10000);
  generateFirstGenomeSegment(g);

  size_t generated_reads = 0;
  size_t current_position = 0;
  size_t real_position = 0;
  size_t remaining_genome = g.length;
  
  while (generated_reads < M) {

    size_t d = generateInterReadDistance();

    // this is artificial however for reasonable values of parameters it should
    // never happen otherwise we woul need a different way of online generating
    // the genome.
    // More specifically if that happens it means that 'd' is higher then a whole
    // genome segment (which should be no less than 10000 in practical cases) for
    // reasonable values of N and M this event will have probability zero for
    // all practical situations and artifically skipping over such 'extreme' values
    // of d will not appreciably change final results
    if (d > (g.length - m - 1)) {
      continue;
    }

    // if we do not have enough generated genome for another read we generate
    // a new genome segment
    if (remaining_genome < d + m) {            
      generateNewGenomeSegment(g);
      std::cout << '\t' << current_position;
      current_position = current_position + m - g.length;
      std::cout << '\t' << current_position << '\n';
      remaining_genome = g.length;
      std::cout << "+-+-+-+-\n";
    }

   
    current_position += d;
    remaining_genome -= d;

    Read current = generateOnlineRead(g.genome,current_position);
    std::cout << current.r << '\t' << current_position << '\n';
    generated_reads++;
    real_position += d;
    current.j = real_position;
  }
}


int main(int argc, char** argv) {
  std::cout << std::endl;
  
  // Important NOT invert (init requires argument to be parsed)
  parseArguments(argc,argv);
  initSimulator();

  if (Options::opts.online) {
    onlineSimulation();
  } else {  
    offlineSimulation();
  }
  
  
  std::cout << std::endl;
  return 0;
}
