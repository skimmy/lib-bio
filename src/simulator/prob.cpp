#include "common.h"

#include <iostream>

#include <boost/math/special_functions/binomial.hpp>


// probability of equality of two calls of the SAME reference position
double p_equal_calls = 0.0;

// after initializations this matrix will contain the probability of tw random
// generated reads showing a overlap of s bases with dh discrepancies
double** p_ran_read_olap;

void initProbabilities() {
// convenience variables
  double p_err = Options::opts.pe;
  double q_err = 1.0 - p_err;
  size_t N = Options::opts.N;
  size_t m = Options::opts.m;

  // init p_equal_calls
  p_equal_calls = q_err * q_err + (1.0 / 3.0) * ( p_err * p_err );


  // init rand_read matrix
  p_ran_read_olap = new double*[m];
  for (size_t i = 0;  i < m; ++i) {
    p_ran_read_olap[i] = new double[m+1];
    for (size_t j = 0; j <= i; ++j) {
      double binom = boost::math::binomial_coefficient<double>(i,j);
      double p_s = (pow(p_equal_calls,j) * pow(1.0 - p_equal_calls, i - j)  * binom ) /  ((double)N * pow(4,i) );
      double p_not_s = 0.0;
      p_ran_read_olap[i][j] = p_s + p_not_s;
      std::cout << p_ran_read_olap[i][j] << '\t';
    }
    std::cout << '\n';
  }
}


void clearProbabilities() {
for (size_t i = 0; i < Options::opts.m; ++i) {
delete[] p_ran_read_olap[i];
  }
  delete[] p_ran_read_olap;
}
