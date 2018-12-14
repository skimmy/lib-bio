/**
 * This is a simulator developed ad-hoc to allow future optimization
 * regardless the changes to other components of the library (which
 * are designed ti be part of a library rather then efficient stand
 * alon tools.
 */

#include <include/common.hpp>
#include <include/options.hpp>

// Basic includes
#include <include/generator.hpp>
#include <include/log.hpp>
#include <include/prob.hpp>
#include <include/util.hpp>


// Task includes
#include <include/align.hpp>
#include <include/chain.hpp>
#include <include/edit.hpp>
#include <include/edit_estimates.hpp>
#include <include/online.hpp>
#include <include/lcs.hpp>

// standard includes
#include <cstdlib>
#include <ctime>

#include <unordered_map>
#include <iostream>
#include <fstream>
#include <algorithm>

// TODO Convert to use proper ns
using namespace lbio::sim;
using namespace lbio::sim::generator;
using namespace lbio::sim::edit;
using namespace lbio::sim::lcs;
using namespace lbio::sim::log;

// output quantities (common to online and offline)
double p_fail = 0.0;
size_t holes = 0;
size_t actually_produced_reads = 0;
double scoreSum = 0.0;
EmpiricalDistribution scoreDist(0,1,10);

// output quantities for oracle simulations

/* 
   An estimation point contains for each s=0,...,m
     - sum of all scores
     - sum of all denominators
     - sum of all denominators
     - number of observations of 's'
   This can then be used to compute average quantities and compare
   them (i.e., esitmated) with the theoretical one (possibly
   approximated).
*/ 
struct EstimationPoint
{  
  double sumNum;
  double sumDen;
  double sumScore;
  int count;
  int hammDist;
};
EstimationPoint * oraclePoints;

EditDistanceSimOutput* edOut;

void initSimulator() {
  initUtil(Options::opts.m);
  initRandomGenerator();
  initProbabilities();
  initChainMatrix();
    
  scoreDist
    = EmpiricalDistribution(0,1,Options::opts.empiricalDistributionStep);
  oraclePoints = new EstimationPoint[Options::opts.m + 1];
  
  for (size_t i = 0; i < Options::opts.m; ++i) {
    oraclePoints[i].sumNum = 0.0;
    oraclePoints[i].sumDen = 0.0;
    oraclePoints[i].sumScore = 0.0;
    oraclePoints[i].count = 0;
    oraclePoints[i].hammDist = 0;
  }

  // set the precision from options
  std::cout.precision(Options::opts.floatPrecision);

  edOut = new EditDistanceSimOutput();
}

void clearSimulator() {
  delete edOut;
  delete[] oraclePoints;
  clearChainMatrix();
  clearProbabilities();
  clearUtil();
}

void outputResults() {
  if (Options::opts.task == Task::EditDist) {
    // CDF ouput requested (-D <distfile> option)
    if (!Options::opts.outputDistribution.empty() &&
	Options::opts.subTask != EDIT_DISTANCE_SUBTASK_SCRIPT_DIST) {
      if (edOut->distPDF) {
	size_t n = Options::opts.N;	
	std::ofstream ofs(Options::opts.outputDistribution, std::ofstream::out);
	logInfo("Writing Edit Distance distribution on "
		+ Options::opts.outputDistribution);
	for (size_t i = 0; i <= n; ++i) {
	  ofs << edOut->distPDF[i] << std::endl;
	}
	ofs.close();
      }
    }
    return;
  }
  if (Options::opts.pipeline) {
    if (Options::opts.task == Task::Oracle) {
      double appNumDen[2];
      for (size_t i = 0; i < Options::opts.m+1; ++i) {
	approximatedScore(i, appNumDen);
	std::cout << i << "\t" << oraclePoints[i].sumScore << "\t"
		  << oraclePoints[i].count << "\t" << oraclePoints[i].sumNum
		  << "\t" << oraclePoints[i].sumDen << "\t" <<appNumDen[0]
		  << "\t" << appNumDen[1] << "\t" << oraclePoints[i].hammDist
		  << "\n";
      }
      return;
    }
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

void recordScoreWithOverlap(double sc, size_t s) {
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
  
  ref = new char[N];
  generateIIDGenome(N,ref);
  std::string s(ref);

  // priority queue is used with position as key so that while
  // extractin reads at once (i.e., emptying the queue) reads will be
  // presented in ordered by position on the reference sequence
  std::priority_queue<Read> reads;
  generateOfflineReads(s, reads, m, Options::opts.M, Options::opts.pe);

  // Temporary variables to count the number of holes, in the future a
  // more sophisticated way (e.g., finite state machine) should be
  // used.
  bool onHole = false;
  
  Read r1 = reads.top();
  reads.pop();
  while(!reads.empty()) {
    Read r2 = reads.top();
    reads.pop();
    size_t s = m - (r2.j - r1.j);
    evaluateChainRelation(r1, r2, s);

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
    }
    r1 = r2;
  }
  
  delete[] ref;
}

void onlineSimulation() {

  bool onHole = false;

  size_t N = Options::opts.N;
  size_t m = Options::opts.m;
  double pe = Options::opts.pe;
  
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

    // this is artificial however for reasonable values of parameters
    // it should never happen otherwise we woul need a different way
    // of online generating the genome.  More specifically if that
    // happens it means that 'd' is higher then a whole genome segment
    // (which should be no less than 10000 in practical cases) for
    // reasonable values of N and M this event will have probability
    // zero for all practical situations and artifically skipping over
    // such 'extreme' values of d will not appreciably change final
    // results
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
   
    Read current = generateOnlineRead(g.genome,current_position, m, pe);
    actual_M++;
    current.j = real_position;

    // here the probabilities are computed and accumulated
    if (prev_read.j != (size_t)-1) {       
      
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

void oracleSimulation() {  
  size_t n = 2 * Options::opts.m;
  size_t m = Options::opts.m;
  double pe = Options::opts.pe;
  double alpha = 1.0 / ((double)Options::opts.N - 2.0 * m + 1);
  double numDen[2];
  
  char* genome = new char[n];  
  // Oracle simulation loops to produce exactly M-1 consecutive pairs
  for (size_t i = 0; i < Options::opts.M - 1; ++i) {
    // generate 2m bases of genome
    generateIIDGenome(n, genome);

    // generate first reads at position 0
    Read r1 = generateOnlineRead(genome, 0, m, pe);
    // generate inter-arrival d
    size_t d = generateInterReadDistance();

    // generate second reads at position d
    if (d >= m) {      
      oraclePoints[0].sumNum += 1.0;
      oraclePoints[0].sumDen += (1.0 / alpha);
      oraclePoints[0].sumScore += alpha;
      oraclePoints[0].count++;
    } else {
      Read r2 = generateOnlineRead(genome, d, m, pe);
      size_t s = m - d;
      oraclePoints[s].hammDist += prefixSuffixHammingDistance(r1.r, r2.r, s);
      double sc = scoreExt(r1.r, r2.r, s,numDen);
      oraclePoints[s].sumNum = numDen[1];
      oraclePoints[s].sumDen = numDen[0];
      oraclePoints[s].sumScore += sc;
      oraclePoints[s].count++;
      
    }
  }
  delete[] genome;
}

void
editDistanceOpMode() {
    

  using AlgorithmBand  = EditDistanceBandApproxLinSpace<lbio_size_t, std::string>;
  using AlgorithmExact = EditDistanceWF<lbio_size_t, std::string>; 
  // The default edit distance mode is
  // Sample
  // Linear Alg
  // mean and variance output
  // no script
  // no sample matrix
  int flags = Options::opts.optFlags; // -f
  size_t n = Options::opts.N;
  int task = Options::opts.subTask; // -B
  edOut->distPDF = new double[n+1];
  std::fill_n(edOut->distPDF, n+1, 0);

  // TASK - Scripts Generation (8)
  if (task == EDIT_DISTANCE_SUBTASK_SCRIPT_DIST) { // -B 8
    logInfo("Task 'Script Distribution'");

    AlgorithmExact alg {n, n, {1,1,1}};
    std::vector<std::string> allScripts {};    
    generate_scripts(n, n, Options::opts.k, allScripts, alg);

    // if file is given save there otherwise use std out
    if (!Options::opts.outputDistribution.empty()) {
      std::ofstream ofs(Options::opts.outputDistribution, std::ofstream::out);
      for (std::string script : allScripts) {
	ofs << script << "\n";
      }
      ofs.close();
    }
    else {
      for (std::string script : allScripts) {
	std::cout << script << "\n";
      }
    }     
    return;
  } 

  // TASK - Algorithms comparison (32)
  if (task == EDIT_DISTANCE_SUBTASK_COMPARE_ALGS) {
    logInfo("Task 'Algorithms comparison'");
    compare_edit_distance_algorithms(n, n, Options::opts.k);
    return;
  }

  // TASK - Edit distance calculation (0, default)
  if (flags & EDIT_DISTANCE_DIFF_BOUNDED_ERROR) { // -f 64
    logInfo("Task 'g(n) Esitmation'");
    size_t k_max = Options::opts.k;
    double precision = Options::opts.precision;
    double z_confidence = Options::opts.confidence;
    // lambda for output (opt_str for extra algorithm output)
    auto print_cb =
      [=](const SampleEstimates& est_n, const SampleEstimates& est_n_2,
	  std::string opt_str = "") {
      std::cout << (n>>1) << "\t" << est_n_2 << "\n"      
                << n << "\t" << est_n  << "\t" << opt_str <<"\n"; };

    lbio_size_t T = static_cast<lbio_size_t>(std::floor(n / 2.0));
    lbio_size_t Tmin = static_cast<lbio_size_t>(std::sqrt(n));
    // Approximation is required find 'optimal' T >= sqrt(n)
    if (Options::opts.approxLevel == 1) { // -A 1 --> Band with optimal width estimation
      if (flags & EDIT_DISTANCE_BANDWIDTH_ESTIMATE) {	
	logInfo("Estimation of optimal bandwidth...");
	Tmin = std::max(1, Options::opts.approxLevel);
	T = optimal_bandwidth(n, precision / 2, 16, Tmin);
	logDebug(debug_string("[Band Estimate]\t",  "~T*: " + std::to_string(T)));
      } else {
	T = Options::opts.approxLevel;
      }      
      AlgorithmBand alg(n, n, T, {1,1,1});
      logInfo("Estimation...");
      std::vector<SampleEstimates> est =
        edit::difference_estimate(n, precision, z_confidence,
				  k_max, alg, print_cb );
    } // -A 1

    if (Options::opts.approxLevel == 2) { // -A 2 Adaptive T
      logInfo("Adaptive Estimation...");
      std::vector<SampleEstimates> est =
	edit::difference_estimate_adaptive(n, precision, z_confidence,
					 k_max, print_cb );     
    }
    return;
  }

  if (flags & EDIT_DISTANCE_BOUNDED_ERROR) { // -f 32
    double precision = Options::opts.precision;
    double z_confidence = Options::opts.confidence;
    SampleEstimates beEst =
      editDistanceErrorBoundedEstimates(n, precision, z_confidence);
    std::cout << beEst.sampleSize << "\t" << beEst.sampleMean << "\t"
	      << beEst.sampleVariance << "\n";
    return;
  }
  // Task Exhaustive 
  if (flags & EDIT_DISTANCE_ESTIMATE_EXHAUSTIVE) {
    if (flags & EDIT_DISTANCE_INFO_PARTIAL) { // -f 5
      logWarning("only \033[1;37mqudratic algorithm\033[0m" 
		 " available with exhaustive option");
      logInfo("Exhaustive edit distance with min-max info");
      edit_distance_exhastive_with_info(n);
    }
    else {
      // Exhasutve (only quadratic)
      logWarning("only \033[1;37mqudratic algorithm\033[0m" 
		 " available with exhaustive option");
      logInfo("Exhaustive edit distance");
      double avgDist = test_exhaustive_edit_distance_encoded(n, edOut->distPDF);
      std::cout << avgDist << std::endl;
    }
  }
  
  else { // Not exhaustive
    // SAMPLE
    size_t k = Options::opts.k;

    if (flags & EDIT_DISTANCE_ALGORITHM_QUADRATIC) {      
      // QUADRATIC + Sample
      if (flags & EDIT_DISTANCE_INFO_PARTIAL) { 
	// Sample + Quadratic + Partial Info
	if (flags & EDIT_DISTANCE_INFO_SCRIPT) { // -f 14
	  logInfo("Sample quadratic algorithms info");
	  AlgorithmExact algExact {n, n, {1,1,1}};
	  EditDistanceSample<AlgorithmExact> generator {n, n};
	  for (size_t i = 0; i < k; ++i) {
	    generator(algExact);
	    auto info = algExact.backtrack();
	    std::cout << info.n_sub << "\t" << info.n_ins
		      << "\t" << info.n_del << "\n";
	  }
	}
	else {
	  logWarning("Sample quadratic info without script not available");
	}
      }
    }

    
    else { 
      // LINEAR + Sample     
      
      if (flags & EDIT_DISTANCE_INFO_PARTIAL) { // -f 4
	// PARTIAL INFO + Sample + Linear
	logWarning("Partial info for linear under developement");
	std::unique_ptr<EditDistanceInfo[]> samples =
	  editDistSamplesInfoLinSpace(n,k);
	
	// auto -> std::unique_ptr<double[]>
	auto subSamples = extractSubstitutionArray(samples.get(), k);
	auto delSamples = extractDeletionArray(samples.get(), k);
	auto insSamples = extractInsertionArray(samples.get(), k);
	
	SampleEstimates subEst
	  = estimatesFromSamples<double>(subSamples.get(), k);
	SampleEstimates delEst
	  = estimatesFromSamples<double>(delSamples.get(), k);
	SampleEstimates insEst
	  = estimatesFromSamples<double>(insSamples.get(), k);

	// If 'verbose' is set  all samples are printed
	if (Options::opts.verbose) {
	  for (size_t i = 0; i < k; ++i) {
	    std::cout << samples[i] << std::endl;
	  }
	}
	std::cout << subEst.sampleMean << "\t" << subEst.sampleVariance << "\n";
	std::cout << delEst.sampleMean << "\t" << delEst.sampleVariance << "\n";
	std::cout << insEst.sampleMean << "\t" << insEst.sampleVariance << "\n";	
      }

      if (flags & EDIT_DISTANCE_MINMAX_STRING) { // -f 16
	logInfo("Sampling periodic/constant string distribution");

	std::cout << std::endl;
	// Sampling e("ACGT^*, Y)
	std::string s(n, 'N');
	for (lbio_size_t i = 0; i<n; ++i) {
	  s[i] = bases[i%4];
	}
	std::vector<size_t> v = edit_samples_fixed_string(n, k, s);
	auto est = estimatesFromSamples(v.cbegin(), v.cend(), k);
	std::cout << "ACGT*\t" << est << "\n";

	// sampling  e("A*", Y)
	s = std::string(n, 'A');
	v = edit_samples_fixed_string(n, k, s);
	est = estimatesFromSamples(v.cbegin(), v.cend(), k);
	std::cout << "A*\t" << est << "\n";
	std::cout << std::endl;	
	return;
      }
      
      else { // -f 0
	// MINIMAL INFO (mean + var) + Sample + Linear
	logInfo("Basic sampling");
	AlgorithmBand alg(n, n, std::ceil(n/2.0), {1,1,1});
	auto samples = edit::edit_distance_samples(n,k, alg);
	auto estimators = estimatesFromSamples<size_t>(samples.get(), k);
	std::cout << estimators.sampleMean << std::endl;
	std::cout << estimators.sampleVariance << std::endl;
      }
    }
  }
}

template <typename _IterT>
lbio_size_t
hamming_distance(_IterT s1_b, _IterT s1_e, _IterT s2_b) {
  lbio_size_t hd = 0;
  for (; s1_b != s1_e; ++s1_b, ++s2_b) {
    hd += static_cast<lbio_size_t>(*s1_b != *s2_b);
  }
  return hd;
}

template <typename _IterT>
lbio_size_t
shift_distance(_IterT s1_b, _IterT s1_e, _IterT s2_b, _IterT s2_e, lbio_size_t n) {
  lbio_size_t sd = hamming_distance(s1_b, s1_e, s2_b);
  lbio_size_t n_2 = n >> 1;
  _IterT s1_m = s1_b;
  _IterT s2_m = s2_b;
  for (lbio_size_t t = 1; t <=n_2; ++t) {
    s1_m++;
    s2_m++;
    lbio_size_t sl = 2*t + hamming_distance(s1_m, s1_e, s2_b);
    lbio_size_t sr = 2*t + hamming_distance(s2_m, s2_e, s1_b);

    sd = std::min(sd, std::min(sl, sr));
  }
  return sd;
}

const lbio_size_t n0_ = 64;
lbio_size_t small_dp_m[n0_][n0_];

lbio_size_t
small_ed(const std::string& s1, const std::string& s2,
	 lbio_size_t i1, lbio_size_t i2, lbio_size_t n) {
  for (lbio_size_t i = 0; i <= n; ++i) {
    small_dp_m[0][i] = i;
    small_dp_m[i][0] = i;
  }
  for (lbio_size_t i = 1; i <= n; ++i) {
    for (lbio_size_t j = 1; j <= n; ++j) {
      lbio_size_t delta = (s1[i1+i-1] == s2[i2+j-1]) ? 0 : 1;
      small_dp_m[i][j] = std::min(std::min(small_dp_m[i-1][j] + 1, small_dp_m[i][j-1] + 1),
				  small_dp_m[i-1][j-1] + delta);
    }
  }
  return small_dp_m[n][n];
}

lbio_size_t
shift_distance_recursive(const std::string& s1, const std::string& s2,
			 lbio_size_t i1, lbio_size_t i2, lbio_size_t n) {
  /*if (n == 1) {
    return static_cast<lbio_size_t>(s1[i1] != s2[i2]);
    }*/
  if (n <= 8) {
    return small_ed(s1, s2, i1, i2, n);
  }
  lbio_size_t sd = shift_distance(s1.begin()+i1, s1.begin()+i1+n, s2.begin()+i2, s2.begin()+i2+n, n);
  lbio_size_t srec_1 = shift_distance_recursive(s1, s2, i1, i2, n >> 1);
  lbio_size_t srec_2 = shift_distance_recursive(s1, s2, i1 + (n >> 1), i2 + (n >> 1), n >> 1);
  return std::min(sd, srec_1 + srec_2);
}

void shift_distance_prototyping() {
   for(lbio_size_t l = 4; l < 17; ++l) {
    lbio_size_t n = 1 << l;
      lbio_size_t k = 200; 
      std::string s1(n, 'N');
      std::string s2(n, 'N');
      lbio_size_t sum_dist = 0;
      for (lbio_size_t i = 0; i < k; ++i) {
	generator::generateIIDString(s1);
	generator::generateIIDString(s2);
	lbio_size_t d = shift_distance_recursive(s1, s2, 0, 0, n);	
	sum_dist += d;
      }
      std::cout << n << "\t" << ((double)sum_dist) / ((double)k*n) << "\n";
  }  
}

//--------------------------------------------------------------------------------

template<typename V_>
void
decrement(V_& v, lbio_size_t k, const V_& def) {
  lbio_size_t i = k-1;
  while(v[i] == 0 and i > 0) {
    v[i] = def[i];
    i--;
  }
  if (i>0) {
    v[i]--;
  }
}

template<typename V_>
void
increment(V_& v, const V_& def, lbio_size_t n_v) {
  lbio_size_t carry = 0;
  lbio_size_t i = n_v-1;
  do {
    v[i] = (v[i]+1) % def[i];
    carry = (v[i] != 0) ? 0 : 1;
    i--;
  } while (carry > 0 && i>=0);
}

template<typename V_>
void
minus_one_non_negative(V_& v, lbio_size_t v_n) {
  for (lbio_size_t i = 0; i < v_n; ++i) {
    v[i] = std::max<lbio_size_t>(0, v[i]-1);
  }
}

void
prototype_partitions() {
    lbio_size_t n = 4;
  std::vector<lbio_size_t> f(n+1,0);
  double N = std::pow(4,n);
  double eed = exhaustive_edit_distance_improved(n, f);
  std::cout << "alpha(" << n << ") = " << eed << "\n";
  for (lbio_size_t d = 0; d <= n; ++d) {
    std::cout << d << "\t" << f[d] << "\t" << (double)f[d]/N << "\n";
  }
  std::cout << "\n\n";
  std::string Sigma = "0123";
  std::vector<lbio_size_t> part;
  part.push_back(2);
  part.push_back(1);
  part.push_back(1);
  part.push_back(1);
  minus_one_non_negative(part, part.size());
  std::vector<lbio_size_t> def(part.size(), 0); 
  
  lbio_size_t total_strings = 1;
  lbio_size_t str_len = 0;
  for (lbio_size_t s = 0; s<Sigma.size(); ++s) {
    def[s] = std::pow(s+1, part[s]);
    total_strings *= def[s];
    str_len += part[s];
  }

  for (lbio_size_t i = 0; i<part.size(); ++i) {
    std::cout << part[i] << "\t" << def[i] << "\n";
  }
  std::cout << total_strings << "\t" << str_len << "\n";

  std::string X(str_len, 'N');
  for (lbio_size_t i = 0; i < part.size(); ++i) {
    std::cout << X << "\n";
  }
}



// custom hash version
//using ColumnStateSet = std::unordered_set<uint32_t, hash_32_bit, equal_high_bits>;
template <typename _K, typename _V>
class ColumnStateSpaceT {  
  typedef typename std::unordered_map<_K, _V>   ColumnStateSet;
  typedef typename ColumnStateSet::iterator     iterator;
  typedef typename std::pair<_K, _V>            pair_type;
  
public:
    
private:
  lbio_size_t _n;  
  ColumnStateSet _set;  
  
public:		    
  ColumnStateSpaceT(lbio_size_t n) : _n {n}, _set {} { }

  // if the key exist add mult to the already existing value
  void insert(_K key, _V mult) {
    auto key_pos = _set.find(key);
    if (key_pos != _set.end()) {
      key_pos->second += mult;
    }
    else {
      _set.insert({key, mult});
    }
  }

  void erase(_K key) {
    _set.erase(key);
  }


  lbio_size_t size() {
    return _set.size();
  }

  iterator begin() {
    return _set.begin();
  }

  iterator end() {
    return _set.end();
  }

  
};

using ColumnStateSpace = ColumnStateSpaceT<uint32_t, uint32_t>;

const uint32_t MinusOne = 0x0;
const uint32_t Zero = 0x1;
const uint32_t One = 0x2;

uint32_t constant_column(lbio_size_t n, uint32_t value) {
  uint32_t col = value & 0x3;
  for (lbio_size_t i = 1; i < n; ++i) {
    col <<= 2;
    col += (0x3 & value);    
  }
  return col;
}

std::string state_to_string(uint32_t state, lbio_size_t n) {
  std::string out {""};
  std::string decoded[] = {"-1", "0", "1", "*"};
  for (lbio_size_t i = 0; i < n; ++i) {
    out = decoded[(state & 0x3)] + ", " + out;
    state >>= 2;
  }
  return "[" + out + "]";
}

void
print_state(ColumnStateSpace& _set, lbio_size_t n) {
  for (auto it_ = _set.begin(); it_ != _set.end(); it_++) {
    std::cout << "(" << state_to_string(it_->first, n) << " [0x" << std::hex << it_->first << "], " <<  std::dec << it_->second << ")\n";
  }
}



// TODO: currently the Hamming mask is a double pointer matrix -> change
ColumnStateSpace
refresh_space(ColumnStateSpace& old_state, lbio_size_t j, lbio_size_t n,
	      lbio_size_t ** h_mask , lbio_size_t sigma) {
  ColumnStateSpace new_state(n);
  auto it_ = old_state.begin();
  lbio_size_t a,b,c,d;
  while(it_ != old_state.end()) {
    // for each mask...
    for (lbio_size_t s_ = 0; s_ < sigma; ++s_) {
      uint32_t Mj = (*it_).first;
      uint32_t Mj1 = 0x0;
      // first element
      a = j;
      b = j+1;
      c = a + (Mj>>2*(n-1) & 0x3) - 1;
      d = std::min(std::min(b+1, c+1), a + h_mask[0][s_]);
      Mj1 += 0x3 & (d-b+1);
      for (lbio_size_t i = 1; i < n; ++i) {
	Mj <<= 2;
	Mj1 <<= 2;
	a = c;
	b = d;
	c = a + (Mj>>2*(n-1-j) & 0x3) - 1;
	d = std::min(std::min(b+1, c+1), a + h_mask[i][s_]);
	Mj1 += 0x3 & (d-b+1);
      }
      new_state.insert(Mj1, (*it_).second);
    }
    it_++;
  }
  return new_state;
}

void
create_hamming_mask(const std::string& str_, const std::string& alphabet_, lbio_size_t** mask) {
  lbio_size_t sigma = alphabet_.size();
  lbio_size_t n = str_.size();
  for (lbio_size_t s = 0; s < sigma; ++s) {
    for (lbio_size_t i = 0; i < n; ++i) {
      mask[i][s] = str_[i] == alphabet_[s] ? 0 : 1;
    }
  }
}



void
prototyping() {
  // !!! WARNING: possibly don't remove next two lines !!!
  std::string proto_task_msg = make_bold("Improved exhaustive");
  logInfo("Working on " + proto_task_msg + " prototyping");
  // -------------------------------------------------------
  //prototype_partitions();
  
  lbio_size_t n = 8;
  lbio_size_t ** mask = allocMatrix<lbio_size_t>(n, 4);
  create_hamming_mask("GACGTAAT", "ACGT", mask);
  //printMatrix(n,4,mask);
  
  ColumnStateSpace state_space(n);
  state_space.insert(constant_column(n, One), 1);
  std::cout << "j = 0\t\t" << state_space.size() << "\n"; 
  print_state(state_space, n);
  std::cout << "\n\n";
  for (lbio_size_t j = 1; j <= n; ++j) {
    state_space = refresh_space(state_space, j, n, mask, 4);
    std::cout << "j = " << j << "\t\t" << state_space.size() << "\n";
    print_state(state_space, n);
    std::cout << "\n\n";
  }
  freeMatrix(n,4,mask);
}

int main(int argc, char** argv) {   
  // Important NOT invert (init requires argument to be parsed)


#ifdef HAVE_BOOST_PROGRAM_OPTIONS
  parseArgumentsBoost(argc,argv);
#else
  parseArguments(argc,argv);
#endif

  
  initSimulator();

  switch (Options::opts.task) {
  case (Task::Test):
    logWarning("Prototyping, for tests run proper binary");
    prototyping();
    exit(0);
  case (Task::Offline):
    offlineSimulation();
    break;
  case (Task::Online):
    onlineSimulation();
    break;
  case (Task::Oracle):
    oracleSimulation();
    break;
  case (Task::AlignScore):
    evaluateAlignmentScore(Options::opts);
    break;
  case (Task::EditDist): // -O 6
    editDistanceOpMode();
    break;
  default:
    std::cout << "Unrecognized operation mode " <<
      static_cast<int>(Options::opts.task) << "\nAborting..\n";
    exit(1);
  }
  outputResults();
  clearSimulator();

  return 0;
}
