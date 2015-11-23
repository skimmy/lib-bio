#include "common.h"

#include <iostream>

#include <boost/math/special_functions/binomial.hpp>


// probability of equality of two calls of the SAME reference position
double p_equal_calls = 0.0;

// probability of 'sequenced position' which represents the probability of a read
// starting at a given position of the reference
double p_read_start = 0.0;

// after initializations this matrix will contain the probability of two random
// generated reads showing a overlap of s bases with dh discrepancies
double** p_rand_read_olap;

// after initializations this matrix will contain in position [s][d] the probability
// of two consecutive reads showing an overlap of 's' bases with 'd' hamming
// distance
double** p_cons_read_olap;

void initProbabilities() {
// convenience variables
  double p_err = Options::opts.pe;
  double q_err = 1.0 - p_err;
  size_t N = Options::opts.N;
  size_t m = Options::opts.m;

  // init p_equal_calls
  p_equal_calls = q_err * q_err + (1.0 / 3.0) * ( p_err * p_err );

  // init p_read_start
  p_read_start = (double)Options::opts.M / (double)Options::opts.N;

  // init rand_read_olap matrix
  p_rand_read_olap = new double*[m];
  for (size_t i = 0;  i < m; ++i) {
    p_rand_read_olap[i] = new double[m+1];
    for (size_t j = 0; j <= i; ++j) {
      double binom = boost::math::binomial_coefficient<double>(i,j);
      double p_s = (pow(p_equal_calls,j) * pow(1.0 - p_equal_calls, i - j)  * binom ) /  ((double)N * pow(4,i) );
      
      double p_not_s = ((double)N-m+1)* pow((3.0/4.0),j) / ( (double)N * pow(4,i-j) ) ;
      p_rand_read_olap[i][j] = p_s + p_not_s;
    }
  }

  // init cons_read_olap matrix
  p_cons_read_olap  = new double*[m];
  for (size_t i = 0; i < m; ++i) {
    p_cons_read_olap[i] = new double[m+1];
    for (size_t j = 0; j <= i ; ++j) {
      double p_s = 0.0;
      double p_not_s = 0.0;

      p_cons_read_olap[i][j] = p_s + p_not_s;
    }
  }
}


void clearProbabilities() {
for (size_t i = 0; i < Options::opts.m; ++i) {
  delete[] p_rand_read_olap[i];
  delete[] p_cons_read_olap[i];
 }
 delete[] p_rand_read_olap;
 delete[] p_cons_read_olap;
}

/**
 * Computes probability of a random overlap of 's' bases with an hamming distance dh
 */
double randomReadsOverlapProbNoErr(size_t s, size_t dh) {
  double _4s = pow(4,s);
  double ind = dh == 0 ? 1.0 : 0.0;
  double ind_pos = (double)(Options::opts.N - 2 * Options::opts.m + 1);
  return ( (_4s*ind) / (ind_pos + ( _4s*ind ) ) );
}
