#include "common.hpp"

#include <iostream>
#include <random>

#include <boost/math/special_functions/binomial.hpp>

// random number generator object for std random library
std::mt19937 gen;

// probability of equality of two calls of the SAME reference position
double p_equal_calls = 0.0;

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

  // init p_equal_calls
  p_equal_calls = q_err * q_err + (1.0 / 3.0) * ( p_err * p_err );

  // init p_read_start
  p_read_start = (double)Options::opts.M / (double)Options::opts.N;

  // init online inter-read distance distribution (i.e., geomtric)
  geom = std::geometric_distribution<>(p_read_start);
}


void clearProbabilities() {

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
double overlappingStringsSum(const std::string & s1, const std::string& s2) {
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
