#include "common.hpp"

#include "generator.hpp"
#include "options.hpp"
#include "prob.hpp"
#include "util.hpp"
#include "edit.hpp"

#include <iostream>
#include <fstream>
#include <map>
#include <memory>

#include <cstring>
#include <cmath>



// ----------------------------------------------------------------------
//                               INFO CONVERSION
// ----------------------------------------------------------------------

std::unique_ptr<double[]> extractSubstitutionArray(const EditDistanceInfo* v, size_t k) {
  std::unique_ptr<double[]> o(new double[k]);
  for (size_t i = 0; i < k; ++i) {
    o[i] = v[i].n_sub;
  }
  return o;
}


std::unique_ptr<double[]> extractDeletionArray(const EditDistanceInfo* v, size_t k) {
  std::unique_ptr<double[]> o(new double[k]);
  
  for (size_t i = 0; i < k; ++i) {
    o[i] = v[i].n_del;
  }
  return o;
}


std::unique_ptr<double[]> extractInsertionArray(const EditDistanceInfo* v, size_t k)  {
  std::unique_ptr<double[]> o(new double[k]);
  for (size_t i = 0; i < k; ++i) {
    o[i] = v[i].n_ins;
  }
  return o;
}

// ----------------------------------------------------------------------
//                        EDIT DISTANCE COMPUTATION
// ----------------------------------------------------------------------

/**
 * \brief Conputes the edit distance between strings s1 and s2 using only
 * linear space (the vecotors passed as parameters). Vectors must be at
 * least m+1 long where m is the length of the second string s2
 */
size_t
editDistanceLinSpace(const std::string& s1, const std::string& s2, size_t* v0, size_t* v1) {
  size_t n1 = s1.size();
  size_t n2 = s2.size();
  size_t n_max = std::max(n1, n2);
  for (size_t i = 0; i < n_max+1; ++i) {
    v0[i] = i;
  }

  for (size_t i = 1; i <= n1; ++i) {
    v1[0] = i;
    for (size_t j = 1; j <= n2; ++j) {
      size_t delta = (s1[i-1] == s2[j-1]) ? 0 : 1;
      v1[j] = std::min( std::min( v0[j] + 1, v1[j-1] + 1), v0[j-1] + delta );
    }
    size_t * tmp = v0;
    v0 = v1;
    v1 = tmp;
  }
  return v0[n2];
}

EditDistanceInfo editDistanceLinSpaceInfo(const std::string& s1, const std::string& s2, EditDistanceInfo* v0, EditDistanceInfo* v1, EditDistanceInfo** sampleMat) {
  size_t n1 = s1.size();
  size_t n2 = s2.size();
  size_t n_max = std::max(n1, n2);

  EditDistanceInfo* tmp = NULL;

  for (size_t i = 0; i <= n_max; ++i) {
    v0[i].n_ins = i;
    v0[i].n_del = 0;
    v0[i].n_sub = 0;
  }

  for (size_t i = 1; i <= n1; ++i) {
    v1[0].n_ins = 0;
    v1[0].n_del = i;
    v1[0].n_sub = 0;
    for (size_t j = 1; j <= n2; ++j) {
      size_t delta = (s1[i-1] == s2[j-1]) ? 0 : 1;
      size_t a = v0[j-1].distance() + delta; // a = M[i-1][j-1] + delta
      size_t b = v0[j].distance() + 1;       // b = M[i-1][j]   + 1     
      size_t c = v1[j-1].distance() + 1;     // c = M[i][j-1]   + 1

      // Case of match
      if (s1[i-1] == s2[j-1]) {
	v1[j] = v0[j-1];
	//	continue;
      } else {
	// Case of substitution NOT WORSE than others
	if ( (a <= b) && (a <= c)) {
	  v1[j] = v0[j-1];
	  v1[j].n_sub++;	
	} else {
	  // In case of tie (equality) always select a deletion
	  if ( b <= c ) {
	    v1[j] = v0[j];
	    v1[j].n_del++;
	  } else {
	    v1[j] = v1[j-1];
	    v1[j].n_ins++;
	  }
	}
      }
      if (sampleMat != NULL) {
	sampleMat[i-1][j-1] += v1[j];
      }
    }
    
    tmp = v0;
    v0 = v1;
    v1 = tmp;

  }
  //  printVector<EditDistanceInfo>(v0,n2+1,"\t"); std::cout << std::endl;
  return v0[n2];
}

// returns the edit distance between strings encoded in two bits form on the 64
// for bits input integers (strings can't be longer than 32 characters). The
// actual lengths of the strings are given as parameters
size_t
editDistanceEncoded(uint64_t s1, size_t n1, uint64_t s2, size_t n2, size_t** dpMatrix) {

  for (size_t i = 1; i < n1+1; ++i) {
    for(size_t j = 1; j < n2+1; ++j) {
      uint64_t x = ( s1 >> 2*(i-1) ) & 0x3; // pre compute matrix {A,C,G,T} x [1...n]
      uint64_t y = ( s2 >> 2*(j-1) ) & 0x3;
      size_t delta = (x == y) ? 0 : 1; // try to find an alternative not involving if
      
      dpMatrix[i][j] = std::min( std::min(dpMatrix[i-1][j]+1, dpMatrix[i][j-1]+1) , dpMatrix[i-1][j-1] + delta ) ;
    }
  }
  return dpMatrix[n1][n2];
}


void
computeAverageDPMatrix(double** dpMatrix, size_t n, size_t m) {
  for (size_t i = 0; i <= n; ++i) {
    dpMatrix[i][0] = i;
  }
  for (size_t j = 0; j <= m; ++j) {
    dpMatrix[0][j] = j;
  }
  for (size_t i = 1; i <= n; ++i) {
    for (size_t j = 1; j <= m; ++j) {
      double minMatch =    std::min( std::min ( dpMatrix[i-1][j] + 1, dpMatrix[i][j-1] + 1), dpMatrix[i-1][j-1] );
      double minMismatch = std::min( std::min ( dpMatrix[i-1][j] + 1, dpMatrix[i][j-1] + 1), dpMatrix[i-1][j-1] + 1);
      dpMatrix[i][j] = 0.25 * minMatch + 0.75 * minMismatch;
    }
  }
}

void
editInfoCompute(EditDistanceInfo& info) {
  info.n_sub = 0;
  info.n_ins = 0;
  info.n_del = 0;
  for (char c : info.edit_script) {
    if (c == 'I') {
      info.n_ins++;
    }
    if (c == 'D') {
      info.n_del++;      
    }
    if (c == 'S') {
      info.n_sub++;
    }
  }
}

void
editDistanceMat(const std::string& s1, const std::string& s2, size_t** dpMatrix) {
  size_t n = s1.size();
  size_t m = s2.size();

  // initialization of first row and column
  for (size_t i = 0; i < n+1; ++i) {
    dpMatrix[i][0] = i;
  }
  for (size_t j = 0; j < m+1; ++j) {
    dpMatrix[0][j] = j;
  }
  
  for (size_t i = 1; i < n+1; ++i) {
    for(size_t j = 1; j < m+1; ++j) {
      size_t delta = (s1[i-1] == s2[j-1]) ? 0 : 1;
      dpMatrix[i][j] = std::min( std::min(dpMatrix[i-1][j]+1, dpMatrix[i][j-1]+1) , dpMatrix[i-1][j-1] + delta ) ;
    }
  }

}

size_t
editDistance(const std::string& s1, const std::string& s2) {
  size_t n = s1.size();
  size_t m = s2.size();
  size_t** dpMatrix = new size_t*[n+1];
  for (size_t i = 0; i < n+1; ++i) {
    dpMatrix[i] = new size_t[m+1];
  }

  editDistanceMat(s1, s2, dpMatrix);
 
  size_t dist = dpMatrix[n][m];
  for (size_t i = 0; i < n+1; ++i) {
    delete[] dpMatrix[i];   
  }
  delete[] dpMatrix;
  return dist;
}

// ----------------------------------------------------------------------
//                    EDIT DISTANCE APPROXIMATIONS
// ----------------------------------------------------------------------

void
editDistanceBandwiseApproxMat(const std::string& s1, const std::string& s2, size_t T, size_t** dpMatrix) {
  size_t n = s1.size();
  size_t m = s2.size();
  size_t INF = n+m+1;
  // init matrix assuming T < min{n,m} (the strictness is crucial to
  // initialize the 'border' diagonals  
  dpMatrix[0][0] = 0;	   
  for (size_t i = 1; i <= T; ++i) {
    dpMatrix[i][0] = i;
  }
  for (size_t j = 1; j <= T; ++j) {
    dpMatrix[0][j] = j;
  }

  for (size_t t = 0; t <= m - (T+1); ++t) {
    dpMatrix[t][T+1+t] = INF;
  }

  for (size_t t = 0; t <= n - (T+1); ++t) {
    dpMatrix[T+1+t][t] = INF;
  }  

  for (size_t i = 1; i <= n; ++i) {
    size_t j_min = (size_t)std::max<int>(1, (int)(i-T));
    size_t j_max = (size_t)std::min<int>(m, (int)(i+T));
    for (size_t j = j_min; j <= j_max; ++j) {      
      size_t delta = (s1[i-1] == s2[j-1]) ? 0 : 1;
      dpMatrix[i][j] = std::min( dpMatrix[i-1][j-1] + delta,
			    std::min(dpMatrix[i-1][j] + 1, dpMatrix[i][j-1] + 1));
    }
  }
}

size_t
editDistanceBandwiseApprox(const std::string& s1, const std::string& s2, size_t T) {
  size_t n = s1.size();
  size_t m = s2.size();
  size_t** dpMatrix = allocMatrix<size_t>(n+1, m+1);
  editDistanceBandwiseApproxMat(s1, s2, T, dpMatrix);  
  size_t dist = dpMatrix[n][m];
  freeMatrix<size_t>(n+1, m+1, dpMatrix);
  return dist;
}


// ----------------------------------------------------------------------
//                        EDIT DISTANCE SAMPLING
// ----------------------------------------------------------------------


void
editDistanceEstimations(size_t n_min, size_t n_max, size_t n_step, size_t k_max) {
  
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
}


std::unique_ptr<size_t[]>
editDistSamples(size_t n, size_t k_samples) {
  std::unique_ptr<size_t[]> v(new size_t[k_samples]);
  std::string s1(n,'N');
  std::string s2(n,'N');
  size_t* v0 = new size_t[n+1];
  size_t* v1 = new size_t[n+1];
  for (size_t k = 0; k < k_samples; ++k) {
    generateIIDString(s1);
    generateIIDString(s2);
    v[k] = editDistanceLinSpace(s1,s2,v0,v1);
  }
  delete[] v0;
  delete[] v1;
  return v;
}

std::unique_ptr<EditDistanceInfo[]>
editDistSamplesInfo(size_t n, size_t k_samples) {
  std::unique_ptr<EditDistanceInfo[]> infos(new EditDistanceInfo[k_samples]);
  std::string s1(n,'N');
  std::string s2(n,'N');

  size_t** dpMatrix = allocMatrix<size_t>(n,n);

  for (size_t k = 0; k < k_samples; ++k) {
    generateIIDString(s1);
    generateIIDString(s2);
    editDistanceMat(s1, s2, dpMatrix);
    editDistanceBacktrack(dpMatrix, s1, s2, infos[k]);
    editInfoCompute(infos[k]);
  }
  
  freeMatrix<size_t>(n, n, dpMatrix);
  
  return infos;
}

std::unique_ptr<EditDistanceInfo[]>
editDistSamplesInfoLinSpace(size_t n, size_t k_samples,  EditDistanceInfo** sampleMat) {
  std::unique_ptr<EditDistanceInfo[]> samples(new EditDistanceInfo[k_samples]);
  std::string s1(n, 'N');
  std::string s2(n, 'N');

  EditDistanceInfo* v0 = new EditDistanceInfo[n+1];
  EditDistanceInfo* v1 = new EditDistanceInfo[n+1];

  for (size_t k = 0; k < k_samples; ++k) {
    generateIIDString(s1);
    generateIIDString(s2);
    samples[k] = editDistanceLinSpaceInfo(s1,s2, v0, v1, sampleMat);
  }

  delete[] v1;
  delete[] v0;
  
  return samples;
}


size_t
sampleEditDistanceDistribution(size_t n, size_t* v0, size_t* v1) {
  std::string s1(n, 'N'); 
  std::string s2(n, 'N');
  generateIIDString(s1);
  generateIIDString(s2);
  size_t d = editDistanceLinSpace(s1, s2, v0, v1);
  return d;
}


// -----------------------------------------------------------------------------
//                            EXHASUTIVE AND BACKTRACK
// -----------------------------------------------------------------------------


double
testExhaustiveEditDistanceEncoded(size_t n, double* freq) {
  size_t** dpMatrix = new size_t*[n+1];
  for (size_t i = 0; i < n+1; ++i) {
    dpMatrix[i] = new size_t[n+1];
  }


  // initialization of first row and column
  for (size_t i = 0; i < n+1; ++i) {
    dpMatrix[i][0] = i;
  }
  for (size_t j = 0; j < n+1; ++j) {
    dpMatrix[0][j] = j;
  }

  uint64_t N = pow(4,n);
  double ed = 0;
  size_t dist = 0;
  for (uint64_t i = 0; i < N; ++i) {
    freq[0]++;
    for (uint64_t j = i+1; j <N; ++j) {
      dist = editDistanceEncoded(i, n, j, n, dpMatrix);
      freq[dist] += 2.0;
      ed += 2*dist;
    }
  }
  for (size_t i = 0; i < n+1; ++i) {
    delete[] dpMatrix[i];
  }
  delete[] dpMatrix;
  return ((double)ed) / ((double) (N*N));
}


void
editDistanceBacktrack(size_t** dpMatrix,const std::string& s1, const std::string& s2, EditDistanceInfo& info) {
  info.edit_script = "";
  size_t i = s1.size();
  size_t j = s2.size();  
  while( i > 0 && j > 0) {
    size_t a = dpMatrix[i-1][j-1];
    size_t b = dpMatrix[i-1][j];
    size_t c = dpMatrix[i][j-1];
    // Match case
    if (s1[i-1] == s2[j-1]) {
      info.edit_script = "M" + info.edit_script;
      i--; j--;
      continue;
    }
    // Substitution case
    if (a <= b && a <=c) {
      info.edit_script = "S" + info.edit_script;
      i--; j--;
      continue;
    }
    // In case of a tie prefer deletion
    if (b <= c) {
      info.edit_script = "D" + info.edit_script;
      i--;
    } else {
      info.edit_script = "I" + info.edit_script;
      j--;
    }      
    
  }

  while (i > 0) {
    info.edit_script = "D" + info.edit_script;
    i--;
  }

  while(j > 0) {
    info.edit_script = "I" + info.edit_script;
    j--;
  }
}

void
closestToDiagonalBacktrack(size_t n, size_t m, size_t** dpMatrix, EditDistanceInfo& info) {
  info.n_sub = 0;
  info.n_del = 0;
  info.n_ins = 0;
  info.edit_script = "";

  size_t i = n;
  size_t j = m;
  size_t a = 0, b = 0, c = 0, d = 0;
  // backtrack until one edge is reached
  while(i > 0 && j > 0) {
    a = dpMatrix[i-1][j-1];
    b = dpMatrix[i-1][j];
    c = dpMatrix[i][j-1];
    d = dpMatrix[i][j];
    
    // Match
    if ( (a == d) && ( a <= b) && (a <= c)) {
      info.edit_script = "M" + info.edit_script;
      i--; j--;
      continue;
    }
    
    // Substitution
    if ( (a < b) && (a < c) ) {
      info.n_sub++;
      info.edit_script = "S" + info.edit_script;
      i--; j--;
      continue;
    }
    // Deletion
    if ( (b < a) && (b < c) ) {
      info.n_del++;
      info.edit_script = "D" + info.edit_script;
      i--;
      continue;
    }
    // Insertion
    if ( (c < a) && (c < b) ) {
      info.n_ins++;
      info.edit_script = "I" + info.edit_script;
      j--;
      continue;
    }
    // Tie between all operations
    if ( (a == b) && (b == c) ) {
      if (i < j) {
	info.n_ins++;
	info.edit_script = "I" + info.edit_script;
	j--;
	continue;
      }
      if (i == j) {
	info.n_sub++;
	info.edit_script = "S" + info.edit_script;
	i--; j--;
	continue;
      }
      if (i > j) {
	info.n_del++;
	info.edit_script = "D" + info.edit_script;
	i--;
	continue;
      }
    }

    // Tie between sub and del
    if (a == b) {

      if (i > j) {
	info.n_del++;
	info.edit_script = "D" + info.edit_script;
	i--;
      } else {
	info.n_sub++;
	info.edit_script = "S" + info.edit_script;
	i--; j--;
      }
      continue;	
    }
    // Tie between sub and ins
    if (a == c) {
      if (i < j) {
	info.n_ins++;
	info.edit_script = "I" + info.edit_script;
	j--;
      } else {
	info.n_sub++;
	info.edit_script = "S" + info.edit_script;
	i--; j--;
      }
      continue;
    }

    // Tie between del and ins
    if (b == c) {
      if (j <= i) {
	info.n_del++;
	info.edit_script = "D" + info.edit_script;
	i--;
      } else {
	info.n_ins++;
	info.edit_script = "I" + info.edit_script;	
	j--;
      }
      continue;

    }
  }

  // backtrack the left edge (if needed)
  while(i > 0) {
    info.n_del++;
    info.edit_script = "D" + info.edit_script;
    i--;
  }

  // bacltrack the top edge (if needed)
  while(j > 0) {
    info.n_ins++;
    info.edit_script  = "I" + info.edit_script;
    j--;
  }
  

}


void editDistanceWithInfo(const std::string& s1, const std::string& s2, EditDistanceInfo& info) {
  size_t n = s1.size();
  size_t m = s2.size();
  size_t** dpMatrix = new size_t*[n+1];
  for (size_t i = 0; i < n+1; ++i) {
    dpMatrix[i] = new size_t[m+1];
  }

  editDistanceMat(s1, s2, dpMatrix);
  editDistanceBacktrack(dpMatrix, s1, s2, info);
  editInfoCompute(info);

  for (size_t i = 0; i < n+1; ++i) {
    delete[] dpMatrix[i];
  }
  delete[] dpMatrix;
}

// ----------------------------------------------------------------------
//                          ESTIMATION PROCEDURE
// ----------------------------------------------------------------------



SampleEstimates
editDistanceErrorBoundedEstimates(size_t n, double precision, double z_delta) {

  size_t* v0 = new size_t[n+1];
  size_t* v1 = new size_t[n+1];
  
  size_t k = 1;
  size_t k_max = Options::opts.k;
  size_t sample = sampleEditDistanceDistribution(n, v0, v1);
  double mean_k = sample;
  double var_k = 0;

  double cumulative_sum = sample;
  double cumulative_quad_sum = sample * sample;
  
  while(k < k_max) {
    k++;
    sample = sampleEditDistanceDistribution(n, v0, v1);
    cumulative_sum += sample;
    cumulative_quad_sum += (sample * sample);
    mean_k = cumulative_sum / ((double)k);
    var_k = ( cumulative_quad_sum - k*(mean_k*mean_k)  ) / ((double)(k-1));
    if (var_k * ( z_delta*z_delta ) < ((double)k) * ( precision * precision )) {
      break;
    }
  }

  delete[] v0;
  delete[] v1;

  SampleEstimates est;
  est.sampleSize = k;
  est.sampleMean = mean_k;
  est.sampleVariance = var_k;

  return est;
}

SampleEstimates
editDistanceRelativeErrorEstimates(size_t n, double e_model, double precision, double z_delta) {

  size_t* v0 = new size_t[n+1];
  size_t* v1 = new size_t[n+1];
  
  size_t k = 1;
  size_t k_max = Options::opts.k;
  size_t k_min = 16;
  size_t sample = sampleEditDistanceDistribution(n, v0, v1);
  double mean_k = sample;
  double var_k = 0;

  double cumulative_sum = sample;
  double cumulative_quad_sum = sample * sample;

  double rho_k = 0;

  do {
    k++;
    sample = sampleEditDistanceDistribution(n, v0, v1);
    cumulative_sum += sample;
    cumulative_quad_sum += (sample * sample);
    mean_k = cumulative_sum / ((double)k);
    var_k = ( cumulative_quad_sum - k*(mean_k*mean_k)  ) / ((double)(k-1));
    rho_k = std::sqrt( var_k / ((double)k));
  } while( k < k_min || (k < k_max &&
			 ( std::abs(mean_k - e_model) < ( rho_k * z_delta / precision ) ) ) );
   
  delete[] v0;
  delete[] v1;

  SampleEstimates est;
  est.sampleSize = k;
  est.sampleMean = mean_k;
  est.sampleVariance = var_k;

  return est;
}

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

// Computes 2 e(n/2) - e(n) with error < (precison) * (value)
std::vector<SampleEstimates>
differenceBoundedRelativeErrorEstimate(size_t n, double precision, double z_delta, size_t k_max) {

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

// ----------------------------------------------------------------------
//                      SCRIPT DISTRIBUTION FUNCTIONS
// ----------------------------------------------------------------------

void
scriptDistributionMatrix(size_t n, size_t m, size_t k, size_t** distMatrix, std::vector<std::string>* scripts) {
  size_t ** dpMatrix = allocMatrix<size_t>(n+1,m+1);
  for (size_t i = 0; i <= n; ++i) {
    for (size_t j = 0; j <= m; ++j) {
      distMatrix[i][j] = 0;
    }
  }

  std::string s1(n,'N');
  std::string s2(n,'N');
  EditDistanceInfo info;

  for (size_t l = 0; l < k; ++l) {
    generateIIDString(s1);
    generateIIDString(s2);
    editDistanceMat(s1, s2, dpMatrix);
    closestToDiagonalBacktrack(n, m, dpMatrix, info);
    // refresh the distribution matrix
    size_t i = 0, j = 0;
    distMatrix[i][j]++;
    for (size_t t = 0; t < info.edit_script.size(); ++t) {
      char c = info.edit_script[t];
      switch(c) {
      case 'M':
	i++; j++;
	break;
      case 'S':
	i++; j++;
	break;
      case 'I':
	j++;
	break;
      case 'D':
	i++;
	break;
      }
      distMatrix[i][j]++;
    }

    if (scripts != nullptr) {
      scripts->push_back(info.edit_script);
    }

  }
  printMatrix<size_t>(n+1,m+1, distMatrix);
  freeMatrix<size_t>(n+1, m+1, dpMatrix);
}

// ----------------------------------------------------------------------
//                        ALGORITHMS COMPARISON
// ----------------------------------------------------------------------

class AlgorithmComparisonResult {
public:
  
  ~AlgorithmComparisonResult() {
    if (exactAlg) {
      delete exactAlg;
    }    
    if (s1) {
      delete s1;      
    }
    if (s2) {
      delete s2;
    }
  }

  void addExact(const EditDistanceInfo& info) {
    if (this->exactAlg) {
      delete this->exactAlg;
    }
    this->exactAlg = new EditDistanceInfo(info);
  }

  void addBandApprox(const EditDistanceInfo& info, size_t T) {
    this->bandApproxAlg[T] = info;
  }

  bool hasExact() {
    return (this->exactAlg != nullptr);
  }

  bool hasBandApproxWithT(size_t T) {
    return (this->bandApproxAlg.count(T));
  }

  EditDistanceInfo getExact() {
    return *(this->exactAlg);
  }

  EditDistanceInfo getBandApproxWithT(size_t T) {
    return this->bandApproxAlg[T];
  }

private:
  EditDistanceInfo* exactAlg = nullptr;
  std::map<size_t, EditDistanceInfo> bandApproxAlg;

  // if needed we may want to return also the string used to calculate
  // to compute the edit distances. Pointersa re used to minimize the
  // required memory (destructor will take care of freeing memory).
  std::string* s1 = nullptr;
  std::string* s2 = nullptr;

};

void
compareEditDistanceAlgorithms(size_t n, size_t m, size_t k, std::ostream& os) {
  size_t T_max = n / 2;
  size_t T_min = 1;

  GeometricProgression<size_t> geom(2, T_min);
  std::vector<size_t> Ts = geom.valuesLeq(T_max);
  Ts.push_back(0);
  

  
  size_t** dpMatrix = allocMatrix<size_t>(n+1, m+1);
  std::vector< std::shared_ptr<AlgorithmComparisonResult> > results;
  std::string s1(n, 'N');
  std::string s2(m, 'N');
  for (size_t l = 0; l < k; ++l) {
    std::shared_ptr<AlgorithmComparisonResult> res =
      std::make_shared<AlgorithmComparisonResult>();
    generateIIDString(s1);
    generateIIDString(s2);

    EditDistanceInfo tmp;
    editDistanceMat(s1, s2, dpMatrix);
    closestToDiagonalBacktrack(s1.size(), s2.size(), dpMatrix, tmp);
    res->addExact(tmp);

    // Approximation for all values of T
    for (size_t T : Ts) {
      editDistanceBandwiseApproxMat(s1, s2, T, dpMatrix);
      closestToDiagonalBacktrack(s1.size(), s2.size(), dpMatrix, tmp);
      res->addBandApprox(tmp, T);
    }
    results.push_back(res);
  }

  os << n << "\t";
  for (size_t T : Ts) {
    os << T << "\t";
  }
  os << std::endl;
  for (auto pRes : results) {
    os << pRes->getExact() << "\t";
    for (size_t T : Ts) {
      os << pRes->getBandApproxWithT(T) << "\t";
    }
    os << std::endl;
    }
  
  freeMatrix<size_t>(n+1, m+1, dpMatrix);
}

// ----------------------------------------------------------------------
//                        OUTPUT RESULT CLASS
// ----------------------------------------------------------------------


EditDistanceSimOutput::~EditDistanceSimOutput() {
  if (this->distPDF) {
    delete[] this->distPDF;
    this->distPDF = NULL;
  }
}
