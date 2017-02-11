/*! \file edit_estimates.hpp 
 */
#ifndef SIM_EDIT_ESTIMATES_H
#define SIM_EDIT_ESTIMATES_H

#include <include/common.hpp>

#include <include/options.hpp>

#include <include/generator.hpp>
#include <include/prob.hpp>

#include <include/edit.hpp>

#include <cmath>

namespace lbio { namespace sim { namespace edit {

/**
   \brief Computes a matrix containing in (i,j) the number of times a
   'closest diagonal' path traverses the cell (i,j)

   \tparam Algorithm The algorithm used for edit distance (must supply
                     \c backtrack() method returning \c EditDistanceInfo).
  
   \param n     Number of rows of the matrix
   \param m     Number of columns of the matrix
   \param k     Number of samples to be used
   \param a     \c std::vector where scritps will be inserted
   \param alg   The instance of \c Algorithm used for calculation and
                backtrack

   \note The \c Algorithm type is used to instantiate a \c
   EditDistanceSample object first and to calculate the backtrack
   afterward. It is therefore required that \c Algorithm type
   implements \c calculate(...) and \c backtrack()

 */
template <class Algorithm>
void
generate_scripts(lbio_size_t n, lbio_size_t m, lbio_size_t k,
		 std::vector<std::string>& scripts, Algorithm& alg) {
  EditDistanceSample<Algorithm> generator {n, m};
  for (; k>0; --k) {
    generator(alg);
    scripts.push_back(alg.backtrack().edit_script);
  }
}



/*!
   \brief Computes the edit distance for \c k_samples pairs of random
   strings with length n. The distances are stored in an array whose
   pointer is returned.
  
   \tparam Algorithm  Algorithm used to compute edit distance

   \param  n          The length of the strings
   \param  k_samples  The number of pairs to compute
  
   \return a  Pointer to an array containing all the distances
 */
template<class Algorithm>
std::unique_ptr<size_t[]>
edit_distance_samples(size_t n, size_t k_samples, Algorithm& alg) {
  std::unique_ptr<size_t[]> v(new size_t[k_samples]);
  std::string s1(n,'N');
  std::string s2(n,'N');
  for (size_t k = 0; k < k_samples; ++k) {
    generator::generateIIDString(s1);
    generator::generateIIDString(s2);
    v[k] = alg.calculate(s1, s2);
  }
  return v;
}


/*!
  \brief Computes \f$2e(n/2) - e(n)\f$ with error < (precison) * (value)
  
  \tparam Algorithm   Class for edit distance algorithm.
  \tparam Func        Callback function invoked with two \c SampleEstimates
  one for n and one for n/2 estimations.
   
  
  \param n          The size of strings to be generated.
  \param precision  The desired _relative precision_.
  \param z_delta    The number of standard deviation for confidence.
  \param k_max      The maximum number of iterations
  \param alg        An instance of `Algorithm` used to calculate the distance
  \param callback   A callable function accepting two \c SampleEstimates
  
  \return A vector \c v containing \c SampleEstimates for e(n/2)
  (\c v[0] )and e(n) (\c v[1] )
 */
template<class Algorithm, typename Func>
std::vector<SampleEstimates> 
difference_estimate(size_t n, double precision, double z_delta,
		   size_t k_max, Algorithm& alg, Func callback) {
  
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
    callback(est_n.toSampleEstimates(), est_n_2.toSampleEstimates());
    // stopping condition
    // - max iteration number reached or
    // - rho < epsilon * |2e(n/2) - e(n)
    
  } while(k < k_max && ( rho >= precision * diff_n / z_delta ));

  return std::vector<SampleEstimates>
    ({ est_n_2.toSampleEstimates(), est_n.toSampleEstimates()}) ;
}


/**
   \overload

   \brief This a convenience proxy function that invokes its
   counterpart with a dummy callback.
 */
template<class Algorithm>
std::vector<SampleEstimates> 
difference_estimate(size_t n, double precision, double z_delta,
		   size_t k_max, Algorithm& alg) {
  
  auto dummy_cb = [](const SampleEstimates& est_n,
		     const SampleEstimates& est_n_2) {};
  return difference_estimate(n, precision, z_delta, k_max, alg, dummy_cb);
}

template <typename _Func>
std::vector<SampleEstimates>
difference_estimate_adaptive(lbio_size_t n, double precision, double z_delta,
			   lbio_size_t kmax, _Func callback) {
  using BandApprAlg = EditDistanceBandApproxLinSpace<lbio_size_t, std::string>;
   std::vector<SampleEstimates> out;
  // lbio_size_t n_2 = static_cast<lbio_size_t>(std::floor(n/2.0));
  // lbio_size_t T = 1;static_cast<lbio_size_t>(std::floor(std::sqrt(n)));
  // lbio::sim::generator::IidPairGenerator gen_n {n, n};
  // lbio::sim::generator::IidPairGenerator gen_n_2 {n_2, n_2};
  // // these will be used for all calculations, only bandwidth value will change
  // BandApprAlg alg_n {n, n, n_2, {1,1,1}};
  // SamplingEstimationProcess est_n {n};
  // SamplingEstimationProcess est_n_2 {n_2};
  // SamplingEstimationProcess appr_est_n {n};
  // SamplingEstimationProcess appr_est_n_2 {n_2};
  // double avg_approximation_error = 1.0;  
  // while(avg_approximation_error > precision / 2.0) {
  //   for (lbio_size_t k = 0; k < k_step; ++k, --kmax) {
  //     auto pair_n = gen_n();      
  //     auto pair_n_2 = gen_n_2();
  //     alg_n.change_bandwidth(n_2);
  //     est_n.newSample(alg_n.calculate(pair_n.first, pair_n.second));
  //     est_n_2.newSample(alg_n.calculate(pair_n_2.first, pair_n_2.second));
  //     alg_n.change_bandwidth(T);
  //     appr_est_n.newSample(alg_n.calculate(pair_n.first, pair_n.second));
  //     appr_est_n_2.newSample(alg_n.calculate(pair_n_2.first, pair_n_2.second));
  //     std::cout << "  " << appr_est_n.sampleMean() - est_n.sampleMean() << "\n";
  //   }

  //   double appr_rho_n = (appr_est_n.sampleMean() / est_n.sampleMean()) - 1;
  //   double appr_rho_n_2 = (appr_est_n_2.sampleMean() / est_n_2.sampleMean()) - 1;
  //   std::cout << avg_approximation_error << "\n";
  //   avg_approximation_error= std::max(appr_rho_n, appr_rho_n_2);
 
  // }  
  // while (kmax > 0) {        
  //   --kmax;
  // }    
   return out;
}

      
} } } // namespaces

#endif