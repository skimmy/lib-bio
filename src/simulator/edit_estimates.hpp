#ifndef SIM_EDIT_ESTIMATES_H
#define SIM_EDIT_ESTIMATES_H

#include "common.hpp"

#include "generator.hpp"
#include "options.hpp"
#include "prob.hpp"
#include "util.hpp"
#include "edit.hpp"

#include <fstream>
#include <cmath>

class SamplingEstimationProcess {
public:
  SamplingEstimationProcess(size_t n)
    : cumulativeSum(0), cumulativeSumSquare(0), k(0)
  {
    this->n = n;
    this->frequency = new size_t[n+1];
    std::fill_n(this->frequency, n+1, 0);
  }
  ~SamplingEstimationProcess() {
    delete[] frequency;
  }

  void newSample(size_t sample) {
    this->frequency[sample]++;
    this->cumulativeSum += sample;
    this->cumulativeSumSquare += (sample * sample);
    this->k++;
  }

  double sampleMean() const {
    return ( (double)cumulativeSum ) / ( (double)k);
  }

  double sampleVariance() const {
    double sMean = sampleMean();
    double meanTerm = ((double)k) * sMean * sMean;
    return ( ( this->cumulativeSumSquare - meanTerm ) / ((double)k-1));
  }

  size_t medianForSampleDistribution() const {
    return medianFromFrequency<size_t>(frequency, n+1);
  }

  size_t sampleSize() const {
    return k;
  }

  void writeFrequencyOnFile(const std::string& path) {
    std::ofstream os(path, std::ofstream::out);
    writeVectorOnStream<size_t>(frequency, n+1, os);
    os.close();
  }

  SampleEstimates toSampleEstimates() const {
    SampleEstimates est;
    est.sampleSize = k;
    est.sampleMean = sampleMean();
    est.sampleVariance = sampleVariance();
    return est;
  }

private:
  size_t n;
  double cumulativeSum;
  double cumulativeSumSquare;
  size_t k;

  size_t* frequency;
};


template<typename Algorithm>
std::vector<SampleEstimates>
differenceBoundedRelativeErrorEstimate(size_t n, double precision, double z_delta, size_t k_max , Algorithm& alg) {
  
  // Generators
  EditDistanceSample<Algorithm> generator_n(n, n);
  EditDistanceSample<Algorithm> generator_n_2(n >> 1, n >> 1);
  
  
  // TODO Need a way to parametrize the output checkpoints
  GeometricProgression<size_t> power_2(2,8);
  LinearProgression<size_t> step_10(8,16);

    

  lbio_size_t k = 1;
  // create sampling process structures (one for 'n' and one for 'n/2' )
  SamplingEstimationProcess est_n(n);
  SamplingEstimationProcess est_n_2(n >> 1);

  // generate new sample
  lbio_size_t sample_n = generator_n(alg);
  lbio_size_t sample_n_2 = generator_n_2(alg);

  // refresh process and variables
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
    // generate new sample
    sample_n = generator_n(alg);
    sample_n_2 = generator_n_2(alg);
  
    // refresh sample processes
    est_n.newSample(sample_n);
    est_n_2.newSample(sample_n_2);
    // get recomputed params...
    mean_n = est_n.sampleMean();
    mean_n_2 = est_n_2.sampleMean();    
    diff_n = 2 * mean_n_2 - mean_n;
    var_n = est_n.sampleVariance();
    var_n_2 = est_n_2.sampleVariance();
    // ...and error terms
    rho = std::sqrt( (4*var_n_2 + var_n) / ((double)k) );

    // middle checkpoint
    /* if (k == step_10.getCurrent() && Options::opts.verbose) {

      std::cout << est_n.sampleSize()
		<< "\t" << mean_n << "\t" << var_n << "\t" << est_n.medianForSampleDistribution()
		<< "\t" << mean_n_2 << "\t" << var_n_2 << "\t" << est_n_2.medianForSampleDistribution()
		<< std::endl;
      step_10.getNext();
    }
    
    if (k == power_2.getCurrent()) {
      
      // print density to file
      est_n.writeFrequencyOnFile("/tmp/density2_" + std::to_string(k) + ".txt");
      power_2.getNext();
 
      }*/
    
    
    // continue looping as long as
    //      rho > threshold
    // and maximum iterations threshold is not reached
  } while(k < k_max && ( rho > precision * diff_n / z_delta ));
    

  std::vector<SampleEstimates> estimates(2);
  estimates[0] = est_n_2.toSampleEstimates();
  estimates[1] = est_n.toSampleEstimates();

  return estimates;
}

/*

// Computes 2 e(n/2) - e(n) with error < (precison) * (value)
std::vector<SampleEstimates>
differenceBoundedRelativeErrorEstimate__DEPRECATED(size_t n, double precision, double z_delta, size_t k_max) {

  // TODO Need a way to parametrize the output checkpoints
  GeometricProgression<size_t> power_2(2,8);
  LinearProgression<size_t> step_10(8,16);

  size_t* v0 = new size_t[n+1];
  size_t* v1 = new size_t[n+1];

  size_t k = 1;
  // create sampling process structures (one for 'n' and one for 'n/2' )
  SamplingEstimationProcess est_n(n);
  SamplingEstimationProcess est_n_2(n >> 1);

  // generate new sample
  size_t sample_n = sampleEditDistanceDistribution(n, v0, v1);
  size_t sample_n_2 = sampleEditDistanceDistribution(n >> 1, v0, v1);

  // refresh process and variables
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
    // generate new sample
    sample_n = sampleEditDistanceDistribution(n, v0, v1);
    sample_n_2 = sampleEditDistanceDistribution(n >> 1, v0, v1);
    // refresh sample processes
    est_n.newSample(sample_n);
    est_n_2.newSample(sample_n_2);
    // get recomputed params...
    mean_n = est_n.sampleMean();
    mean_n_2 = est_n_2.sampleMean();    
    diff_n = 2 * mean_n_2 - mean_n;
    var_n = est_n.sampleVariance();
    var_n_2 = est_n_2.sampleVariance();
    // ...and error terms
    rho = std::sqrt( (4*var_n_2 + var_n) / ((double)k) );

    // middle checkpoint
    if (k == step_10.getCurrent() && Options::opts.verbose) {

      std::cout << est_n.sampleSize()
		<< "\t" << mean_n << "\t" << var_n << "\t" << est_n.medianForSampleDistribution()
		<< "\t" << mean_n_2 << "\t" << var_n_2 << "\t" << est_n_2.medianForSampleDistribution()
		<< std::endl;
      step_10.getNext();
    }
    
    if (k == power_2.getCurrent()) {
      
      // print density to file
      est_n.writeFrequencyOnFile("/tmp/density2_" + std::to_string(k) + ".txt");
      power_2.getNext();
 
    }
    
    
    // continue looping as long as
    //      rho > threshold
    // and maximum iterations threshold is not reached
  } while(k < k_max && ( rho > precision * diff_n / z_delta ));
    
  delete[] v1;
  delete[] v0;
  std::vector<SampleEstimates> estimates(2);
  estimates[0] = est_n_2.toSampleEstimates();
  estimates[1] = est_n.toSampleEstimates();

  return estimates;
}

*/

#endif
