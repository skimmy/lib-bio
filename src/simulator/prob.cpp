#include "common.hpp"

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

void initProbabilities() {

  std::random_device rd;
  gen = std::mt19937(rd());
  
  // convenience variables
  double p_err = Options::opts.pe;
  double q_err = 1.0 - p_err;

  // init p_equal_calls and q_equal_call
  p_equal_calls = q_err * q_err + (1.0 / 3.0) * ( p_err * p_err );
  q_equal_calls = 1.0 - p_equal_calls;

  // init p_read_start
  p_read_start = (double)Options::opts.M / (double)Options::opts.N;

  // init online inter-read distance distribution (i.e., geomtric)
  geom = std::geometric_distribution<>(p_read_start);
}


void clearProbabilities() {

}

double indicatorErr(const std::string& r1, const std::string& r2, size_t s) {
  if (s > 0) {
    size_t hamm_d = prefixSuffixHammingDistance(r1, r2, s);
    return ( pow(p_equal_calls, s - hamm_d) * pow(q_equal_calls, hamm_d) ) ;
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
    sum += pow(4,s) * ( indicator_ab + indicator_ba );
  }
  size_t indic_m = (prefixSuffixHammingDistance(s2,s1,m) == 0) ? 1 : 0;
  sum += indic_m * pow(4,m);
  return sum;
}

double overlappingStringsSumWithErr(const std::string& s1, const std::string& s2) {
  double sum = 0.0;
  size_t m = Options::opts.m;
  double peq = p_equal_calls;
  double qeq = 1.0 - p_equal_calls;
  for (size_t s = 1; s <= m-1; ++s) {
    // WARNING: If needed here we can use lookup tables to make things faster
    double tab = indicatorErr(s1, s2, s);
    double tba = indicatorErr(s2, s1, s);
    sum += pow(4,s) * (tab + tba);
  }
  sum += pow(4,m) * indicatorErr(s1, s2, m);
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


double score(const std::string& r1, const std::string& r2, size_t s) {
  
  double num = 1.0;
  double den = 1.0;
  // error case
  if (Options::opts.pe > 0) {
    den = (double)Options::opts.N - 2.0 * (double)Options::opts.m + 1.0
      + overlappingStringsSumWithErr(r1,r2);
    
    num = indicatorErr(r1, r2, s) *pow(4,s);
    
  }
  // no error case
  else {
    den = (double)Options::opts.N - 2.0 * (double)Options::opts.m + 1.0
      + overlappingStringsSum(r1,r2);
    num = pow(4,s);
  }
  return num / den;
}
