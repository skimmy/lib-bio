#ifndef SIM_PROB_H
#define SIM_PROB_H

#include <include/util.hpp>

#include <vector>
#include <fstream>
#include <iostream>

// Class representing empirical distribution as histogram array
class EmpiricalDistribution {  
public:
  EmpiricalDistribution(double a, double b, size_t N);

  // given a sample x in [a,b] returns the index where x must be counted
  size_t indexForSample(double x) const;
  // returns the relative frequenxy (i.e., count[i] / total) for the i-th index
  double valueAtIndex(size_t i) const;
  // returns the number of elements at the i-th index
  double countAtIndex(size_t i) const { return f[i]; }

  // returns the number of intervals 'n'
  double getIntervalCount() const { return n; }
  
  // add a new sample x to the distribution
  void addSample(double x);

  // returns the emprical cumulative distribution derived from
  // currently stored samples.
  void getCDF(std::vector<double>& cdf)  const;
    
  
private:
  std::vector<double> f;
  size_t n;
  double xa;
  double xb;
  double step; // defined as (b-a)/n
  double total;
};

struct ScoreSumFreq
{
  double sSum;
  int sFreq;
};

class SampleEstimates
{
public:
  size_t sampleSize;
  
  double sampleMean;
  double sampleVariance;

  SampleEstimates() : sampleSize(0), sampleMean(0), sampleVariance(0) {}

  friend std::ostream& operator<<(std::ostream& os,
				  const SampleEstimates& est) {
    os << est.sampleSize << "\t" << est.sampleMean
       << "\t" << est.sampleVariance;
    return os;
  }
};

/**
 * \brief Computes the estimators in the SampleEstimates class from
 * k samples of type T for which conversion to double must be possible
 */

template<typename T>
SampleEstimates estimatesFromSamples(T samples[], size_t k) {
  SampleEstimates est;
  est.sampleSize = k;

  // estimates the 'sample mean'
  for (size_t i = 0; i < est.sampleSize; ++i) {
    est.sampleMean += (double)samples[i];    
  }
  est.sampleMean /= (double)est.sampleSize;

  // estimates the 'sample variance'
  double tmp = 0;
  for (size_t i = 0; i < est.sampleSize; ++i) {
    tmp = (samples[i] - est.sampleMean);
    est.sampleVariance += tmp*tmp;
  }
  est.sampleVariance /= (double)(est.sampleSize - 1);
  return est;
}

/**
 * \brief Computes the estimators for the given range
 */
template<typename T>
SampleEstimates estimatesFromSamples(T begin, T end, lbio_size_t k) {
  SampleEstimates est;
  est.sampleSize = k;

  // estimates the 'sample mean'
  for (T it = begin; it != end; ++it) {
    est.sampleMean += (double)*it;    
  }
  est.sampleMean /= (double)est.sampleSize;

  // estimates the 'sample variance'
  double tmp = 0;

  for (T it = begin; it != end; ++it) {
    tmp = (*it - est.sampleMean);
    est.sampleVariance += tmp*tmp;
  }
  est.sampleVariance /= (double)(est.sampleSize - 1);
  return est;
  
}

void initProbabilities();
void clearProbabilities();
double randomReadsOverlapProbNoErr(const std::string& s1, const std::string& s2, size_t s);
double overlappingStringsSum(const std::string & s1, const std::string& s2);
double overlappingStringsSumWithErr(const std::string& s1, const std::string& s2);
size_t generateInterReadDistance();

size_t percentileIndex(const std::vector<double>& cdf, double perc);

double approximatedScore(size_t s, double* num_den);
double approximatedScore(size_t s);

double score(const std::string& r1, const std::string& r2, size_t s);
double scoreExt(const std::string& r1, const std::string& r2, size_t s, double* num_den);

/*
   \brief Returns the index of the FIRST elements x such that 
       P[X < x] >= 0.5
 */
template<typename T>
size_t medianFromFrequency(T f[], size_t n) {
  // first compute the sum...
  T sum = 0;
  for (size_t i = 0; i < n; ++i) {
    sum += f[i];
  }
  
  // ...then find when 50% of the mass is covered
  T pSum = 0;
  for (size_t i = 0; i < n; ++i) {
    pSum += f[i];
    if ( ((double)pSum) / ( (double)sum) >= 0.5) {
      return i;
    }    
  }
  return n-1;
}

namespace lbio { namespace prob {

/**
 * \brief This class is used to represent sampling processes.
 * 
 * The class provides facility to compute sample mean and sample
 * variance progressively as samples are supplied.
 */
class SamplingEstimationProcess {
public:
  SamplingEstimationProcess(size_t n_);
  ~SamplingEstimationProcess();

  void 
  newSample(size_t sample);

  double 
  sampleMean() const;

  double 
  sampleVariance() const;

  double 
  standardError() const;
  
  size_t 
  sampleSize() const;

  size_t
  medianForSampleDistribution() const;

  void
  writeFrequencyOnFile(const std::string& path);

  SampleEstimates 
  toSampleEstimates() const;

private:
  size_t n;
  double cumulativeSum;
  double cumulativeSumSquare;
  size_t k;

  size_t* frequency;
};

} } // namespaces

#endif
