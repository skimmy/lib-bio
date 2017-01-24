#ifndef SIM_EDIT_ESTIMATES_H
#define SIM_EDIT_ESTIMATES_H

#include "common.hpp"

#include "prob.hpp"
#include "options.hpp"
#include "edit.hpp"

#include <cmath>

namespace lbio
{
namespace sim
{
namespace edit
{


/*
 * @brief Computes k_the edit distance for k_samples pairs of random
 * strings with length n. The distances are stored in an array whose
 * pointer is returned.
 *
 * @tparam Algorithm  Algorithm used to compute edit distance
 * @param  n          The length of the strings
 * @param  k_samples  The number of pairs to compute
 *
 * @return a  Pointer to an array containing all the distances
 */
template<class Algorithm>
std::unique_ptr<size_t[]>
edit_distance_samples(size_t n, size_t k_samples, Algorithm& alg) {
  std::unique_ptr<size_t[]> v(new size_t[k_samples]);
  std::string s1(n,'N');
  std::string s2(n,'N');
  for (size_t k = 0; k < k_samples; ++k) {
    generateIIDString(s1);
    generateIIDString(s2);
    v[k] = alg.calculate(s1, s2);
  }
  return v;
}


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
 *
 * @return A `vector v` containing `SampleEstimates` for e(n/2)
 *         (`v[0]`)and e(n) (`v[1]`)
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



// Old and note used functions
/**
 * Use a Monte-Carlo sampling technique to estimate the edit distance between
 * random strings of same length.
 *
 * \param n_min The minimum length to be tested
 * \param n_max The maximum lengths to be tested
 * \param n_step The incremente on length 
 * \param k_max The number of samples for each length
 * 
 */
/*void
editDistanc eEstimations(size_t n_min, size_t n_max, size_t n_step, size_t k_max) {  
  std::cout << std::endl;
  for (size_t n = n_min; n <= n_max; n += n_step) {
    std::string s1(n,'N');
    std::string s2(n,'N');
    double AED = 0;
    for (size_t k = 1; k <= k_max; ++k) {
      generateIIDString(s1);
      generateIIDString(s2);
      AED += editDistance(s1,s2);
    }
    std::cout << n << "\t" << ( AED / k_max) << std::endl;
  }
  std::cout << std::endl;
  }      */

  
} // ::edit
} // ::sim
} // ::lbio

#endif
