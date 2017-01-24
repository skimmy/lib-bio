#ifndef SIM_EDIT_ESTIMATES_H
#define SIM_EDIT_ESTIMATES_H

#include "common.hpp"

#include "prob.hpp"
#include "options.hpp"
#include "edit.hpp"

#include <fstream>
#include <cmath>

/**
 * @brief Computes 2 e(n/2) - e(n) with error < (precison) * (value)
 *
 * @tparam Algorithm Class for edit distance algorithm.
 *
 * @param n          The size of strings to be generated.
 * @param precision  The desired _relative precision_.
 * @param z_delta    The number of standard deviation for confidence.
 * @param k_max      The maximum number of iterations
 * @param alg        An instance of `Algorithm` used to calculate the distance
 */
template<class Algorithm>
std::vector<SampleEstimates> 
difference_stimate(size_t n, double precision, double z_delta,
		   size_t k_max, Algorithm& alg) {
  
  // Generators
  EditDistanceSample<Algorithm> generator_n(n, n);
  EditDistanceSample<Algorithm> generator_n_2(n >> 1, n >> 1);
  
  lbio_size_t k = 1;
  // create sampling process structures (one for 'n' and one for 'n/2' )
  SamplingEstimationProcess est_n(n);
  SamplingEstimationProcess est_n_2(n / 2);

  // generate new sample and refresh sampling processes
  lbio_size_t sample_n = generator_n(alg);
  lbio_size_t sample_n_2 = generator_n_2(alg);
  est_n.newSample(sample_n);
  est_n_2.newSample(sample_n_2);
  
  // sample means, variances (0 with 1 sample) and diff 2e(n/2) - e(n)
  double mean_n = est_n.sampleMean();
  double mean_n_2 = est_n_2.sampleMean();
  double diff_n = 2 * mean_n_2 - mean_n;
  double var_n = 0;
  double var_n_2 = 0;
  // error term
  double rho = 0;

  do {
    k++;
    // generate new sample and refresh sampling processes
    sample_n = generator_n(alg);
    sample_n_2 = generator_n_2(alg);
    est_n.newSample(sample_n);
    est_n_2.newSample(sample_n_2);

    // recompute params and error
    mean_n = est_n.sampleMean();
    mean_n_2 = est_n_2.sampleMean();    
    diff_n = 2 * mean_n_2 - mean_n;
    var_n = est_n.sampleVariance();
    var_n_2 = est_n_2.sampleVariance();
    rho = std::sqrt( (4*var_n_2 + var_n) / ((double)k) );

    // stopping condition
    // - max iteration number reached or
    // - rho < epsilon * |2e(n/2) - e(n)
  } while(k < k_max && ( rho >= precision * diff_n / z_delta ));
    

  return std::vector<SampleEstimates>
    ({ est_n_2.toSampleEstimates(), est_n.toSampleEstimates()}) ;
}

#endif
