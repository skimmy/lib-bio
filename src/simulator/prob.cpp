#include <include/common.hpp>

#include <include/options.hpp>

#include <include/prob.hpp>
#include <include/util.hpp>

#include <vector>
#include <iostream>
#include <random>

#include <boost/math/special_functions/binomial.hpp>

// random number generator object for std random library
std::mt19937 gen;

// probability of equality of two calls of the SAME reference position
double p_equal_calls = 0.0;
double q_equal_calls = 0.0;

// probability of 'sequenced position' which represents the probability of a read
// starting at a given position of the reference
double p_read_start = 0.0;

// geometric distribution for inter-read distance
std::geometric_distribution<> geom;

// lookup tables for probability of the form
// p^x * (1-p)^{m-s}
double * power_peq_lookup = NULL; // p^i
double * power_qeq_lookup = NULL; // q^i = (1-p)^i

/* lookup tables for the approximated expected score defined as:
                         4^s * tildeI(s, peq)
     tildeE(A,B,s) =  ----------------------------
                      4^s * tildeI(s, peq) + N - 1
 */
double * approxExpScoreNum = NULL;
double * approxExpScoreDen = NULL;

void initProbabilities() {

  std::random_device rd;
  gen = std::mt19937(rd());
  
  // convenience variables
  double p_err = Options::opts.pe;
  double q_err = 1.0 - p_err;

  // init p_equal_calls and q_equal_call
  p_equal_calls = q_err * q_err + (1.0 / 3.0) * ( p_err * p_err );
  q_equal_calls = (1.0 - p_equal_calls) / 3.0;

  // init p_read_start
  p_read_start = (double)Options::opts.M / (double)Options::opts.N;

  // init online inter-read distance distribution (i.e., geomtric)
  geom = std::geometric_distribution<>(p_read_start);

  // lookup tables
  int m = Options::opts.m;
  int N = Options::opts.N;
  power_peq_lookup = new double[m+1];
  power_qeq_lookup = new double[m+1];
  approxExpScoreNum = new double[m+1];
  approxExpScoreDen = new double[m+1];
  for (int s = 0; s <= m; ++s) {
    power_peq_lookup[s] = pow(p_equal_calls, s);
    power_qeq_lookup[s] = pow(q_equal_calls, s);
    double tmp = p_equal_calls * p_equal_calls + q_equal_calls* ( 1 - p_equal_calls );
    double tildeI = pow(tmp, s );
    approxExpScoreNum[s] = pow(4,s) * tildeI;
    approxExpScoreDen[s] = (pow(4,s) * tildeI ) + N - 1;

  }
}


void clearProbabilities() {
  delete[] approxExpScoreNum;
  delete[] approxExpScoreDen;
  delete[] power_qeq_lookup;
  delete[] power_peq_lookup;

}

double indicatorErr(const std::string& r1, const std::string& r2, size_t s) {
  if (s > 0) {
    size_t hamm_d = prefixSuffixHammingDistance(r1, r2, s);
    return (power_peq_lookup[s - hamm_d] * power_qeq_lookup[hamm_d]);
  } 
  return 1.0;
    
}

double indicatorNoErr(const std::string& r1, const std::string& r2, size_t s) {
  return (prefixSuffixHammingDistance(r1, r2, s) == 0 );
}

double probabilityReads(const std::string& s1, const std::string& s2, size_t d) {
  
  // probability of D = d distance between reads
  double p_s = p_read_start * pow(1.0 - p_read_start, d);
  // no overlap
  if (d >= Options::opts.m) {
    return pow(4,2*Options::opts.m) * p_s;
  }

  size_t dh = prefixSuffixHammingDistance(s1,s2,Options::opts.m - d);
  // probability of dh unequal bases on the common part
  double p_err = pow(p_equal_calls, (Options::opts.m - d - dh)) * pow(1.0 - p_equal_calls, dh) ;
  // probability of indipendent parts
  double p_ind = pow(4,2*d);

  return p_s * p_err * p_ind;
}


/**
 * Compute
 *  \sum_{s}{I(A,B,s)4^s}
 */
double overlappingStringsSum(const std::string& s1, const std::string& s2) {
  size_t m = Options::opts.m;
  double sum = 0.0;
  for (size_t s = 1; s <= m-1; ++s) {
    size_t indicator_ab = (prefixSuffixHammingDistance(s2,s1,s) == 0) ? 1 : 0;
    size_t indicator_ba = (prefixSuffixHammingDistance(s1,s2,s) == 0) ? 1 : 0;
    sum += power4_lookup[s] * ( indicator_ab + indicator_ba );
  }
  size_t indic_m = (prefixSuffixHammingDistance(s2,s1,m) == 0) ? 1 : 0;
  sum += indic_m * power4_lookup[m];
  return sum;
}

double overlappingStringsSumWithErr(const std::string& s1, const std::string& s2) {
  double sum = 0.0;
  size_t m = Options::opts.m;
  for (size_t s = 1; s <= m-1; ++s) {
    // WARNING: If needed here we can use lookup tables to make things faster
    double tab = indicatorErr(s1, s2, s);
    double tba = indicatorErr(s2, s1, s);
    sum += power4_lookup[s] * (tab + tba);
  }
  sum += power4_lookup[m] * indicatorErr(s1, s2, m);
  return sum;
}

/**
 * Computes probability of a random overlap of 's' bases with an hamming distance dh
 */
double randomReadsOverlapProbNoErr(const std::string& s1, const std::string& s2, size_t s) {
  size_t dh_for_s = prefixSuffixHammingDistance(s1,s2,s);
  if (dh_for_s != 0) {
    return 0;
  }
  double olap_p = overlappingStringsSum(s1,s2) + (double)Options::opts.N - 2 * (double)Options::opts.m + 1.0;
  return ( pow(4,s) / ( olap_p) ) ;
}

/**
 * Gets the distance until the next read
 */
size_t generateInterReadDistance() {
  return geom(gen);
}

EmpiricalDistribution::EmpiricalDistribution(double a, double b, size_t N)
  : f(N+1,0.0)
{
  this->xa = a;
  this->xb = b;
  this->n = N+1;
  this->step = (b-a) / (double)N;
  this->total = 0;
}

size_t EmpiricalDistribution::indexForSample(double x) const {
  return floor( (x - xa) / step);
}


double EmpiricalDistribution::valueAtIndex(size_t i) const {
  return (f[i] / (double)total);
}
  
void EmpiricalDistribution::addSample(double x) {
  total++;
  size_t i = indexForSample(x);
  f[i] += 1;
}

void EmpiricalDistribution::getCDF(std::vector<double>& cdf) const {
  if (this->n == 0) {
    return;
  }
  cdf.reserve(this->n);
  cdf[0] = this->f[0];
  for (size_t i = 1; i < this->n; ++i) {
    cdf[i] = cdf[i-1] + f[i];
  }
  // normalize between 0 and 1
  for (size_t i = 0; i < this->n; ++i) {
    cdf[i] /= cdf[this->n - 1];
  }
}

size_t percentileIndex(const std::vector<double>& cdf, double perc) {
  size_t idx = 0;
  size_t n = cdf.size();
  for (size_t i = 0; i < n; ++i) {
    if(cdf[i] >= perc) {
      break;
    }
    ++idx;
  }
  return idx;
}

double approximatedScore(size_t s, double* num_den) {
  num_den[0] = approxExpScoreNum[s];
  num_den[1] = approxExpScoreDen[s];
  return approxExpScoreNum[s] / approxExpScoreDen[s];
}

double approximatedScore(size_t s) {
  return approxExpScoreNum[s] / approxExpScoreDen[s];
}

double score(const std::string& r1, const std::string& r2, size_t s) {

  static double iidTerm = power4_lookup[0] *
    ((double)Options::opts.N - 2.0 * (double)Options::opts.m + 1.0);
  
  double den = iidTerm  + overlappingStringsSumWithErr(r1,r2);
  double num = indicatorErr(r1, r2, s) * power4_lookup[s];

  return num / den;
}

double scoreExt(const std::string& r1, const std::string& r2, size_t s, double* num_den) {

  static double iidTerm = power4_lookup[0] *
    ((double)Options::opts.N - 2.0 * (double)Options::opts.m + 1.0);
  
  num_den[0] = iidTerm  + overlappingStringsSumWithErr(r1,r2);
  num_den[1] = indicatorErr(r1, r2, s) * power4_lookup[s];

  return num_den[1] / num_den[0];
}
