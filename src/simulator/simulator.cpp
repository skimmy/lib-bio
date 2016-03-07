/**
 * This is a simulator developed ad-hoc to allow future optimization regardless
 * the changes to other components of the library (which are designed ti be part
 * of a library rather then efficient stand alon tools.
 */

#include "common.hpp"

#include <cstdlib>
#include <ctime>

#include <iostream>


char bases[] = {'A', 'C', 'G', 'T'};
char revBases[128];
Options Options::opts;

// output quantities (common to online and offline)
double p_fail = 0.0;
size_t holes = 0;
size_t actually_produced_reads = 0;
double scoreSum = 0.0;
EmpiricalDistribution scoreDist(0,1,10);

void initSimulator() {
  revBases['A'] = revBases['a'] = 0;
  revBases['C'] = revBases['c'] = 1;
  revBases['G'] = revBases['g'] = 2;
  revBases['T'] = revBases['t'] = 3;
  scoreDist = EmpiricalDistribution(0,1,Options::opts.empiricalDistributionStep);
  srand(time(NULL));
  initProbabilities();
  initChainMatrix();

  // initFalsePositiveMatrix();

}

void clearSimulator() {
  //clearFalsePositiveMatrix();
  clearChainMatrix();
  clearProbabilities();

}

void outputResults() {
  if (Options::opts.pipeline) {
    std::cout << p_fail << std::endl;
  } else {
    std::cout << "P[Fail]    = " << p_fail << std::endl;
    std::cout << "P[Success] = " << 1.0 - p_fail << std::endl;
    std::cout << "#[Holes]   = " << holes << std::endl;
    std::cout << "#[Reads]   = " << actually_produced_reads << std::endl;
  }
  if (!Options::opts.outputDistribution.empty()) {
    std::ofstream ofs(Options::opts.outputDistribution, std::ofstream::out);
    for (size_t i = 0; i < scoreDist.getIntervalCount(); ++i) {
      ofs << scoreDist.valueAtIndex(i) << '\n';
    }
    ofs.close();
  }

  if (!Options::opts.outputCDF.empty()) {
    std::ofstream ofs(Options::opts.outputCDF, std::ofstream::out);
    std::vector<double> cdf(scoreDist.getIntervalCount());
    scoreDist.getCDF(cdf);
    for (size_t i = 0; i < scoreDist.getIntervalCount(); ++i) {
      ofs << cdf[i] << "\n";
    }
    ofs.close();
    std::cout << scoreDist.valueAtIndex(percentileIndex(cdf,0.001)) << "\n";
  }
  
}

void recordScore(double p_ab) {
  p_fail += 1.0 - p_ab;
  scoreSum += p_ab;
  scoreDist.addSample(p_ab);
}

void offlineSimulation() {

  char* ref = NULL;
  
  size_t N = Options::opts.N;
  size_t m = Options::opts.m;
  size_t M = Options::opts.M;
  size_t Nbar = N - m + 1;
  double pe = Options::opts.pe; 
  
  ref = new char[N];
  generateIIDGenome(N,ref);
  std::string s(ref);

  // priority queue is used with position as key so that while extractin reads at
  // once (i.e., emptying the queue) reads will be presented in ordered by
  // position on the reference sequence
  std::priority_queue<Read> reads;
  generateOfflineReads(s, reads);

  // Temporary variables to count the number of holes, in the future a more
  // sophisticated way (e.g., finite state machine) should be used.
  bool onHole = false;
  
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
      recordScore(p_ab);

    } else {
      if (onHole == false) {
	onHole = true;
	holes++;
      }
      addNonOverlapRecord(r2.j - r1.j - m);
      double x = (double)Options::opts.N - 2.0 * (double)Options::opts.m + 1.0
	+ overlappingStringsSum(r1.r, r2.r);
      recordScore(1.0 / x);
      //      p_fail += 1.0 - ( 1.0 / x);
      //      scoreSum += ( 1.0 / x);
    }
    r1 = r2;
  }
  
  delete[] ref;
}

void onlineSimulation() {

  bool onHole = false;

  size_t N = Options::opts.N;
  size_t m = Options::opts.m;
  size_t M = Options::opts.M;  
  
  GenomeSegment g(N, m, MAX_GENOME_SEGMENT_LENGTH);
  generateFirstGenomeSegment(g);

  size_t generated_reads = 0;
  size_t current_position = 0;
  size_t real_position = 0;
  size_t remaining_genome = g.length;

  size_t actual_M = 0;

  Read prev_read("", -1);
  
  while (real_position < N - m) {

    size_t d = generateInterReadDistance();
    real_position += d;

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

    // in this case we need to generate new genome segment...
    if (remaining_genome < m + d) {
      size_t tmp = current_position + d;
      if (tmp < g.length) {
	generateNewGenomeSegment(g, g.length - tmp);
      } else {
	generateNewGenomeSegment(g, 0);
      }
      current_position = 0;
      remaining_genome = g.length;
    } else {
      // ...otherwise we simply update counters
      current_position += d;    
      remaining_genome -= d;
    }
   
    Read current = generateOnlineRead(g.genome,current_position);
    actual_M++;
    current.j = real_position;

    // here the probabilities are computed and accumulated
    if (prev_read.j != -1) {       
      
      if (d > m) {
	if (!onHole) {
	  holes++;
	}
	onHole = true;
	// non-overlap case...
	if (Options::opts.approxLevel < 0) {
	  double sc = score(prev_read.r, current.r , 0);	  
	  recordScore(sc);
	}
		
      } else {
	onHole = false;
	// overlap case...
	size_t s = m - d;
	double sc = score(prev_read.r, current.r , s);
	recordScore(sc);
      }
    }
    
    generated_reads++;
    current.j = real_position;
    prev_read = current;
  }
  actually_produced_reads = actual_M;
}


int main(int argc, char** argv) {   
  // Important NOT invert (init requires argument to be parsed)
  parseArguments(argc,argv);
  initSimulator();

  if (Options::opts.test) {
    testAll();
  } else {

    if (Options::opts.online) {
      onlineSimulation();
    } else {  
      offlineSimulation();
    }
    outputResults();
    
  }
  return 0;
}
