#include "common.hpp"

#include <iostream>
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
//                             DP MATRIX UTILS
// ----------------------------------------------------------------------


size_t**
allocDPMatrix(size_t n, size_t m) {
  size_t** dpMatrix = new size_t*[n+1];
  for (int i = 0; i <= n; ++i) {
    dpMatrix[i] = new size_t[m+1];
  }
  return dpMatrix;
}

void
freeDPMatrix(size_t** dpMatrix, size_t n, size_t m) {
  for (int i = 0; i <= n; ++i) {
    delete[] dpMatrix[i];
  }
  delete[] dpMatrix;
}

void
printDPMatrix(size_t** dpMatrix, size_t n, size_t m) {
  for (size_t i = 0; i <= n; ++i) {
    for (size_t j = 0; j <=m; ++j) {
      std::cout << dpMatrix[i][j] << '\t';
    }
    std::cout << std::endl;
  }
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
  size_t n_max = MAX(n1, n2);
  for (int i = 0; i < n_max+1; ++i) {
    v0[i] = i;
  }

  for (size_t i = 1; i <= n1; ++i) {
    v1[0] = i;
    for (size_t j = 1; j <= n2; ++j) {
      size_t delta = (s1[i-1] == s2[j-1]) ? 0 : 1;
      v1[j] = MIN( MIN( v0[j] + 1, v1[j-1] + 1), v0[j-1] + delta );
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
  size_t n_max = MAX(n1, n2);

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

  for (int i = 1; i < n1+1; ++i) {
    for(int j = 1; j < n2+1; ++j) {
      uint64_t x = ( s1 >> 2*(i-1) ) & 0x3; // pre compute matrix {A,C,G,T} x [1...n]
      uint64_t y = ( s2 >> 2*(j-1) ) & 0x3;
      size_t delta = (x == y) ? 0 : 1; // try to find an alternative not involving if
      
      dpMatrix[i][j] = MIN( MIN(dpMatrix[i-1][j]+1, dpMatrix[i][j-1]+1) , dpMatrix[i-1][j-1] + delta ) ;
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
      double minMatch =    MIN( MIN ( dpMatrix[i-1][j] + 1, dpMatrix[i][j-1] + 1), dpMatrix[i-1][j-1] );
      double minMismatch = MIN( MIN ( dpMatrix[i-1][j] + 1, dpMatrix[i][j-1] + 1), dpMatrix[i-1][j-1] + 1);
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
  
  for (int i = 1; i < n+1; ++i) {
    for(int j = 1; j < m+1; ++j) {
      size_t delta = (s1[i-1] == s2[j-1]) ? 0 : 1;
      dpMatrix[i][j] = MIN( MIN(dpMatrix[i-1][j]+1, dpMatrix[i][j-1]+1) , dpMatrix[i-1][j-1] + delta ) ;
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
  for (int i = 0; i < n+1; ++i) {
    delete[] dpMatrix[i];   
  }
  delete[] dpMatrix;
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
  size_t* v0 = new size_t[n];
  size_t* v1 = new size_t[n];
  for (int k = 0; k < k_samples; ++k) {
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

  size_t** dpMatrix = allocDPMatrix(n,n);

  for (size_t k = 0; k < k_samples; ++k) {
    generateIIDString(s1);
    generateIIDString(s2);
    editDistanceMat(s1, s2, dpMatrix);
    editDistanceBacktrack(dpMatrix, s1, s2, infos[k]);
    editInfoCompute(infos[k]);
  }
  
  freeDPMatrix(dpMatrix, n, n);
  
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


// -----------------------------------------------------------------------------
//                            EXHASUTIVE AND BACKTRACK
// -----------------------------------------------------------------------------


double
testExhaustiveEditDistanceEncoded(size_t n) {
  size_t** dpMatrix = new size_t*[n+1];
  for (int i = 0; i < n+1; ++i) {
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
  for (uint64_t i = 0; i < N; ++i) {
    for (uint64_t j = i+1; j <N; ++j) {
      ed += 2*editDistanceEncoded(i, n, j, n, dpMatrix);      
    }
  }
  for (int i = 0; i < n+1; ++i) {
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
    size_t x = dpMatrix[i][j];
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


void editDistanceWithInfo(const std::string& s1, const std::string& s2, EditDistanceInfo& info) {
  size_t n = s1.size();
  size_t m = s2.size();
  size_t** dpMatrix = new size_t*[n+1];
  for (int i = 0; i < n+1; ++i) {
    dpMatrix[i] = new size_t[m+1];
  }

  editDistanceMat(s1, s2, dpMatrix);
  editDistanceBacktrack(dpMatrix, s1, s2, info);
  editInfoCompute(info);

  for (int i = 0; i < n+1; ++i) {
    delete[] dpMatrix[i];
  }
  delete[] dpMatrix;
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
  double var_k = 0, var_k_1 = 0;

  double cumulative_sum = sample;
  double cumulative_quad_sum = sample * sample;
  
  while(k < k_max) {
    k++;
    sample = sampleEditDistanceDistribution(n, v0, v1);
    cumulative_sum += sample;
    cumulative_quad_sum += (sample * sample);
    mean_k = cumulative_sum / ((double)k);
    var_k_1 = var_k;
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

// Computes 2 e(n/2) - e(n) with error < (precison) * (value)
SampleEstimates
differenceBoundedRelativeErrorEstimate(size_t n, double precision, double z_delta, size_t k_max) {

  size_t* v0 = new size_t[n+1];
  size_t* v1 = new size_t[n+1];

  size_t k = 1;

  size_t sample_n = sampleEditDistanceDistribution(n, v0, v1);
  size_t sample_n_2 = sampleEditDistanceDistribution(n >> 1, v0, v1);

  double cumul_sum_n = sample_n;
  double cumul_sum_n_2 = sample_n_2;

  double cumul_sum_square_n = sample_n * sample_n;
  double cumul_sum_square_n_2 = sample_n_2 * sample_n_2;

  double mean_n = cumul_sum_n / ((double) k);
  double mean_n_2 = cumul_sum_n_2 / ( (double) k);

  double diff_n = 2 * mean_n_2 - mean_n;

  double var_n = ( cumul_sum_square_n - k * mean_n * mean_n ) / ( (double) k-1);
  double var_n_2 = ( cumul_sum_square_n_2 - k * mean_n_2 * mean_n_2) / ( (double)k-1 );

  double rho = 0;

  do {
    k++;

    sample_n = sampleEditDistanceDistribution(n, v0, v1);
    sample_n_2 = sampleEditDistanceDistribution(n >> 1, v0, v1);

    cumul_sum_n += sample_n;
    cumul_sum_n_2 += sample_n_2;

    cumul_sum_square_n += sample_n * sample_n;
    cumul_sum_square_n_2 += sample_n_2 * sample_n_2;

    mean_n = cumul_sum_n / ((double) k);
    mean_n_2 = cumul_sum_n_2 / ( (double) k);

    diff_n = 2 * mean_n_2 - mean_n;

    var_n = ( cumul_sum_square_n - k * mean_n * mean_n ) / ( (double) k-1);
    var_n_2 = ( cumul_sum_square_n_2 - k * mean_n_2 * mean_n_2) / ( (double)k-1 );

    rho = std::sqrt( (4*var_n_2 + var_n) / ((double)k) );

    
  } while(k < k_max && ( rho > precision * diff_n / z_delta ));
    
  delete[] v1;
  delete[] v0;
    
  SampleEstimates est;
  est.sampleMean = diff_n;
  est.sampleVariance = rho;
  est.sampleSize = k;
  return est;
}
