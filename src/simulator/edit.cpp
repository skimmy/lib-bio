#include <include/common.hpp>
#include <include/edit.hpp>

#include <include/options.hpp>

#include <include/prob.hpp>
#include <include/generator.hpp>
#include <include/util.hpp>


#include <fstream>
#include <map>
#include <unordered_map>
#include <future>
#include <thread>
#include <cmath>


using namespace lbio::sim::generator;


//////////////////////////////////////////////////////////////////////
//                    EDIT INFO IMPLE AND HELPERS
//////////////////////////////////////////////////////////////////////

std::unique_ptr<double[]>
extractSubstitutionArray(const EditDistanceInfo* v, size_t k) {
  std::unique_ptr<double[]> o(new double[k]);
  for (size_t i = 0; i < k; ++i) {
    o[i] = v[i].n_sub;
  }
  return o;
}


std::unique_ptr<double[]>
extractDeletionArray(const EditDistanceInfo* v, size_t k) {
  std::unique_ptr<double[]> o(new double[k]);
  
  for (size_t i = 0; i < k; ++i) {
    o[i] = v[i].n_del;
  }
  return o;
}


std::unique_ptr<double[]>
extractInsertionArray(const EditDistanceInfo* v, size_t k)  {
  std::unique_ptr<double[]> o(new double[k]);
  for (size_t i = 0; i < k; ++i) {
    o[i] = v[i].n_ins;
  }
  return o;
}

EditDistanceInfo::EditDistanceInfo()
  : n_sub {0}, n_del {0}, n_ins {0}, edit_script {""}
{ }

EditDistanceInfo::EditDistanceInfo(lbio_size_t s, lbio_size_t d, lbio_size_t i)
  : n_sub {s}, n_del {d}, n_ins {i}, edit_script {""}
{ }

EditDistanceInfo::EditDistanceInfo(lbio_size_t c)
  : EditDistanceInfo(c,c,c)
{ } 

size_t
EditDistanceInfo::distance() const {
  return n_sub + n_ins + n_del;
}
  
void
EditDistanceInfo::reset() {
  n_sub = 0; n_del = 0; n_ins = 0;
}

bool
EditDistanceInfo::operator<(const EditDistanceInfo& i) const {
  return (this->distance() < i.distance());
}
  
bool
EditDistanceInfo::operator==(const EditDistanceInfo& i) const {
  return (n_sub == i.n_sub && n_del == i.n_del && n_ins == i.n_ins);
}
  
bool
EditDistanceInfo::operator!=(const EditDistanceInfo& i) const {
  return !(*this == i);
}
  
EditDistanceInfo&
EditDistanceInfo::operator+=(const EditDistanceInfo& rhs) {
  n_sub += rhs.n_sub;
  n_del += rhs.n_del;
  n_ins += rhs.n_ins;
  return *this;
}

EditDistanceInfo
operator+(EditDistanceInfo lhs, const EditDistanceInfo& rhs) {
  lhs += rhs;
  return lhs;
}

EditDistanceInfo&
EditDistanceInfo::operator*=(lbio_size_t scalar) {
  n_sub *= scalar;
  n_del *= scalar;
  n_ins *= scalar;
  return *this;
}

EditDistanceInfo
operator*(EditDistanceInfo lhs, lbio_size_t rhs) {
  lhs *= rhs;
  return lhs;
}

 
std::ostream&
operator<<(std::ostream& out, const EditDistanceInfo& info) {
  out << info.n_sub << " " << info.n_del << " " << info.n_ins;
  return out;
}  

//////////////////////////////////////////////////////////////////////
//                      EDIT DISTANCE COMPUTATION
//////////////////////////////////////////////////////////////////////

/**
 * \brief Conputes the edit distance between strings s1 and s2 using only
 * linear space (the vecotors passed as parameters). Vectors must be at
 * least m+1 long where m is the length of the second string s2
 */
size_t
editDistanceLinSpace(const std::string& s1, const std::string& s2,
		     size_t* v0, size_t* v1) {
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
      v1[j] = std::min(std::min(v0[j]+1, v1[j-1]+1), v0[j-1]+delta);
    }
    size_t * tmp = v0;
    v0 = v1;
    v1 = tmp;
  }
  return v0[n2];
}

EditDistanceInfo
editDistanceLinSpaceInfo(const std::string& s1, const std::string& s2,
			 EditDistanceInfo* v0, EditDistanceInfo* v1,
			 EditDistanceInfo** sampleMat) {
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
  return v0[n2];
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

//////////////////////////////////////////////////////////////////////
//                    EDIT DISTANCE APPROXIMATIONS
//////////////////////////////////////////////////////////////////////

void
editDistanceBandwiseApproxMat(const std::string& s1, const std::string& s2,
			      size_t T, size_t** dpMatrix) {
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
				 std::min(dpMatrix[i-1][j] + 1,
					  dpMatrix[i][j-1] + 1));
    }
  }
}

size_t
editDistanceBandwiseApprox(const std::string& s1, const std::string& s2,
			   size_t T) {
  size_t n = s1.size();
  size_t m = s2.size();
  size_t** dpMatrix = allocMatrix<size_t>(n+1, m+1);
  editDistanceBandwiseApproxMat(s1, s2, T, dpMatrix);  
  size_t dist = dpMatrix[n][m];
  freeMatrix<size_t>(n+1, m+1, dpMatrix);
  return dist;
}


//////////////////////////////////////////////////////////////////////
//                       EDIT DISTANCE SAMPLING
//////////////////////////////////////////////////////////////////////


std::unique_ptr<EditDistanceInfo[]>
editDistSamplesInfoLinSpace(size_t n, size_t k_samples, std::string alphabet,
			    EditDistanceInfo** sampleMat) {
  std::unique_ptr<EditDistanceInfo[]> samples(new EditDistanceInfo[k_samples]);
  std::string s1(n, 'N');
  std::string s2(n, 'N');

  EditDistanceInfo* v0 = new EditDistanceInfo[n+1];
  EditDistanceInfo* v1 = new EditDistanceInfo[n+1];

  for (size_t k = 0; k < k_samples; ++k) {
    generateIIDString(s1, alphabet);
    generateIIDString(s2, alphabet);
    samples[k] = editDistanceLinSpaceInfo(s1,s2, v0, v1, sampleMat);
  }

  delete[] v1;
  delete[] v0;
  
  return samples;
}

std::vector<size_t>
edit_samples_fixed_string(size_t n, size_t k_samples, const std::string& s2, std::string alphabet) {
  std::vector<size_t> samples;
  std::string s1(n, 'N');

  size_t* v0 = new size_t[n+1];
  size_t* v1 = new size_t[n+1];

  for (size_t k = 0; k < k_samples; ++k) {
    generateIIDString(s1, alphabet);
    samples.push_back(editDistanceLinSpace(s1,s2, v0, v1));
  }

  delete[] v1;
  delete[] v0;
  
  return samples;
}

// ----------------------------------------------------------------------
//                          ESTIMATION PROCEDURE
// ----------------------------------------------------------------------

SampleEstimates
editDistanceErrorBoundedEstimates(size_t n, std::string alphabet, double precision,
				  double z_delta, size_t k_min) {
  using BandApprox = EditDistanceBandApproxLinSpace<lbio_size_t, std::string>;
  BandApprox alg(n, n, std::floor(n/2.0), {1,1,1});
  EditDistanceSample<BandApprox> gen(n, n, alphabet);
  
  size_t k = 1;
  size_t k_max = Options::opts.k;
  size_t sample = gen(alg); 
  double mean_k = sample;
  double var_k = 0;

  double cumulative_sum = sample;
  double cumulative_quad_sum = sample * sample;
  
  while(k < k_max) {
    k++;
    sample = gen(alg);
    cumulative_sum += sample;
    cumulative_quad_sum += (sample * sample);
    mean_k = cumulative_sum / ((double)k);
    var_k = ( cumulative_quad_sum - k*(mean_k*mean_k)  ) / ((double)(k-1));
    // if  ( (var_k * ( z_delta*z_delta ) < ((double)k) * ( precision * precision ))
    // 	  && (k > k_min) ){
    //   break;
    // }
    if  (var_k * ( z_delta*z_delta ) < ((double)k) * ( precision * precision )) {
      break;
    }
  }

  SampleEstimates est;
  est.sampleSize = k;
  est.sampleMean = mean_k;
  est.sampleVariance = var_k;

  return est;
}

SampleEstimates
editDistanceRelativeErrorEstimates(size_t n, std::string alphabet, double e_model,
				   double precision, double z_delta) {

  using BandApprox = EditDistanceBandApproxLinSpace<lbio_size_t, std::string>;
  BandApprox alg(n, n, std::floor(n/2.0), {1,1,1});
  EditDistanceSample<BandApprox> gen(n, n, alphabet);

  
  size_t k = 1;
  size_t k_max = Options::opts.k;
  size_t k_min = 16;
  size_t sample = gen(alg); 
  double mean_k = sample;
  double var_k = 0;

  double cumulative_sum = sample;
  double cumulative_quad_sum = sample * sample;

  double rho_k = 0;

  do {
    k++;
    sample = gen(alg); 
    cumulative_sum += sample;
    cumulative_quad_sum += (sample * sample);
    mean_k = cumulative_sum / ((double)k);
    var_k = ( cumulative_quad_sum - k*(mean_k*mean_k)  ) / ((double)(k-1));
    rho_k = std::sqrt( var_k / ((double)k));
  } while( k < k_min || (k < k_max
			 && ( std::abs(mean_k - e_model)
			      < ( rho_k * z_delta / precision ) ) ) );
   

  SampleEstimates est;
  est.sampleSize = k;
  est.sampleMean = mean_k;
  est.sampleVariance = var_k;

  return est;
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


// --------------------------------------------------------------------
// NAMESPACES BEGIN

namespace lbio { namespace sim { namespace edit {


// returns the edit distance between strings encoded in two bits form
// on the 64 for bits input integers (strings can't be longer than 32
// characters). The actual lengths of the strings are given as
// parameters
lbio_size_t
edit_distance_encoded(uint64_t s1, lbio_size_t n1, uint64_t s2,
		    lbio_size_t n2, lbio_size_t** dpMatrix) {
  for (lbio_size_t i = 1; i < n1+1; ++i) {
    for(lbio_size_t j = 1; j < n2+1; ++j) {
      // pre compute matrix {A,C,G,T} x [1...n]
      uint64_t x = ( s1 >> 2*(i-1) ) & 0x3; 
      uint64_t y = ( s2 >> 2*(j-1) ) & 0x3;
      lbio_size_t delta = (x == y) ? 0 : 1;       
      dpMatrix[i][j] =
	std::min( std::min(dpMatrix[i-1][j]+1, dpMatrix[i][j-1]+1),
		  dpMatrix[i-1][j-1] + delta ) ;
    }
  }
  return dpMatrix[n1][n2];
}

void
edit_distance_exhastive_with_info(lbio_size_t n) {
  lbio_size_t** dpMatrix = allocMatrix<lbio_size_t>(n+1,n+1);
  // initialization of first row and column
  for (lbio_size_t i = 0; i < n+1; ++i) {
    dpMatrix[i][0] = i;
  }
  for (lbio_size_t j = 0; j < n+1; ++j) {
    dpMatrix[0][j] = j;
  }

  uint64_t N = pow(4,n);
  double tot_dist = 0;
  
  double min_dist = n;
  uint64_t min_center = 0;
  double max_dist = 0;
  uint64_t max_center = 0;
  
  for (uint64_t i = 0; i < N; ++i) {
    lbio_size_t dist = 0;
    
    for (uint64_t j = 0; j <N; ++j) {
      dist += edit_distance_encoded(i, n, j, n, dpMatrix);
    }

    double ddist = dist / static_cast<double>(N);
    tot_dist += dist;
    if (ddist < min_dist) {
      min_center = i;
      min_dist = ddist;      
    }
    if (ddist > max_dist) {
      max_center = i;
      max_dist = ddist;
    }
  }
  tot_dist /= static_cast<double>(N*N);
  freeMatrix(n+1, n+1, dpMatrix);
  std::cout << "Min: " << min_center << "\t" << min_dist << "\n";
  std::cout << "Max: " << max_center << "\t" << max_dist << "\n";
  std::cout << "Avg: " << tot_dist << "\n";
}

double
test_exhaustive_edit_distance_encoded(lbio_size_t n, double* freq) {

  lbio_size_t** dpMatrix = allocMatrix<lbio_size_t>(n+1,n+1);

  // initialization of first row and column
  for (lbio_size_t i = 0; i < n+1; ++i) {
    dpMatrix[i][0] = i;
  }
  for (lbio_size_t j = 0; j < n+1; ++j) {
    dpMatrix[0][j] = j;
  }

  uint64_t N = pow(4,n);
  double ed = 0;
  lbio_size_t dist = 0;
  for (uint64_t i = 0; i < N; ++i) {
    freq[0]++;
    for (uint64_t j = i+1; j <N; ++j) {
      dist = edit_distance_encoded(i, n, j, n, dpMatrix);
      freq[dist] += 2.0;
      ed += 2*dist;
    }
  }
  freeMatrix(n+1, n+1, dpMatrix);
  return ((double)ed) / ((double) (N*N));
}

////////////////////////////////////////////////////////////////////////////////
//                          COLUMN STATE EXHAUSTIVE
////////////////////////////////////////////////////////////////////////////////

uint64_t constant_column(lbio_size_t n, uint32_t value) {
  uint64_t col = value & 0x3;
  for (lbio_size_t i = 1; i < n; ++i) {
    col <<= 2;
    col += (0x3 & value);    
  }
  return col;
}

std::string state_to_string(uint64_t state, lbio_size_t n) {
  std::string out {""};
  std::string decoded[] = {"-1", "0", "1", "*"};
  for (lbio_size_t i = 0; i < n; ++i) {
    out = decoded[(state & 0x3)] + ", " + out;
    state >>= 2;
  }
  return "[" + out + "]";
}

std::vector<lbio_size_t> state_to_column(uint64_t state, lbio_size_t j, lbio_size_t n) {
  std::vector<lbio_size_t> column(n+1);
  column[0] = j;
  int vals[] {-1, 0, 1, 0xFF};
  for (lbio_size_t i = 1; i <= n; ++i) {
    lbio_size_t idx = (state >> 2*(n-i)) & 0x3;
    column[i] = column[i-1] + vals[idx];
  }
  return column;
}

void
print_state(ColumnStateSpace& _set, lbio_size_t n) {
  for (auto it_ = _set.begin(); it_ != _set.end(); it_++) {
    std::cout << "(" << state_to_string(it_->first, n) << " [0x" << std::hex
	      << it_->first << "], " <<  std::dec << it_->second << ")\n";
  }
}


ColumnStateSpace
refresh_space(ColumnStateSpace& old_state, lbio_size_t j, lbio_size_t n,
	      lbio_size_t ** h_mask , lbio_size_t sigma) {
  ColumnStateSpace new_state(n);
  auto it_ = old_state.begin();
  lbio_size_t a,b,c,d;
  while(it_ != old_state.end()) {
    // for each mask...
    for (lbio_size_t s_ = 0; s_ < sigma; ++s_) {
      uint64_t Mj = (*it_).first;
      uint64_t Mj1 = 0x0;
      // first element
      a = j;
      b = j+1;
      c = a + (Mj>>2*(n-1) & 0x3) - 1;
      d = std::min(std::min(b+1, c+1), a + h_mask[0][s_]);
      Mj1 += 0x3 & (d-b+1);
      for (lbio_size_t i = 1; i < n; ++i) {
	Mj1 <<= 2;
	a = c;
	b = d;
	c = a + (Mj>>2*(n-1-i) & 0x3) - 1;
	d = std::min(std::min(b+1, c+1), a + h_mask[i][s_]);
	Mj1 += 0x3 & (d-b+1);
      }
      new_state.insert(Mj1, (*it_).second);
    }
    it_++;
  }
  return new_state;
}


ColumnStateSpace
state_space_compute(std::string x, std::string alphabet) {
  lbio_size_t sigma = alphabet.size();
  lbio_size_t ** mask = allocMatrix<lbio_size_t>(x.size(), sigma);
  create_hamming_mask(x, alphabet, mask);
  ColumnStateSpace states(x.size());
  states.insert(constant_column(x.size(), One), 1);
  for (lbio_size_t j = 0; j < x.size(); ++j) {
    states = refresh_space(states, j, x.size(), mask, sigma);
  }
  freeMatrix(x.size(), sigma, mask);  
  return states;
}

template <typename _Iter, typename _Call>
void
state_space_compute_iter(_Iter beg_, _Iter end_, _Call action, std::string alphabet) {
  while(beg_ != end_) {
    ColumnStateSpace states = state_space_compute(*beg_, alphabet);
    action(*beg_, states);
    ++beg_;
  }
}

std::vector<lbio_size_t>
state_space_to_distribution(ColumnStateSpace& states) {
  lbio_size_t n = states.get_n();
  std::vector<lbio_size_t> v(n+1, 0);
  for (auto it_ = states.begin(); it_ != states.end(); it_++) {    
    v[state_to_column(it_->first, n, n)[n]] += it_->second;
  }
  return v;
}

std::string
distribution_to_string(const std::vector<lbio_size_t> v) {
  std::string s {""};
  for(auto it = v.begin(); it != v.end(); it++) {
    s += std::to_string(*it) + ", ";
  }
  return s;
}

std::string get_permutation_invariant_form(std::string s, std::string alphabet) {
  std::string out(s.size(), 'N');
  lbio_size_t i, m = 0;
  std::unordered_map<char, char> map_;
  while(i < s.size()) {
    if (map_.count(s[i]) <= 0) {
      map_[s[i]] = alphabet[m++];
    }
    out[i] = map_[s[i]];
    i += 1;
  }
  return out;
}

template <typename It_>
std::vector<std::pair<std::string, lbio_size_t>>
remove_symmetry_invariant(It_ begin, It_ end, std::string alphabet) {
  std::map<std::string, lbio_size_t> pair_map(begin, end);
  for (; begin != end; ++begin) {
    std::string s = begin->first;
    std::string rev(s.rbegin(), s.rend());
    // CASE 1 - Palindrome strings: nothing todo
    if (s == rev) {
      continue;
    }    
    std::string p_rev = get_permutation_invariant_form(rev, alphabet);
    // CASE 2 - s and its invariant reverse are equal: same as palindrome
    if (s == p_rev) {
      continue;
    }
    // CASE 3 - s and its invariant reverse are different: keep only
    // one and sum multiplicities
    auto itm_ = pair_map.find(p_rev);
    if (itm_ != pair_map.end()) {
      itm_->second += begin->second;
      pair_map.erase(s);
    }
  }
  return std::vector<std::pair<std::string, lbio_size_t>>(pair_map.begin(), pair_map.end());
} 

void
edit_distance_eccentricity(lbio_size_t n, std::ostream& os, std::string alphabet) {
  // iterator for all the strings
  AlphabetIterator it(n, alphabet);
  // for each string compute eccentricoty and write distribution on 'os'
  state_space_compute_iter(it, it.end(), [&os] (std::string x, ColumnStateSpace& space) {
      os << x << ", " << distribution_to_string(state_space_to_distribution(space)) << "\n";
					 }, alphabet);
}

double
eccentricity_for_string(std::string x, std::string alphabet) {
  ColumnStateSpace states = state_space_compute(x, alphabet);
  double sum = 0;
  auto v = state_space_to_distribution(states);
  for (lbio_size_t i = 1; i < v.size(); ++i) {
    sum += i*v[i];
  }  
  return sum/(std::pow(4,x.size())*x.size());
}

template <typename _It>
double
eccentricity_with_symmetries_iterator(_It begin, _It end, lbio_size_t n, std::string alphabet) {
  double sum = 0;
  while(begin != end) {
    ColumnStateSpace states = state_space_compute(begin->first, alphabet);
    auto dist = state_space_to_distribution(states);
    double dist_sum = 0;
    for (lbio_size_t i = 1; i < dist.size(); ++i) {
      dist_sum += i*dist[i];
    }  
    double ecc = dist_sum / std::pow(alphabet.size(), n);
    sum += ecc*begin->second;
    ++begin;
  }
  return sum;
}

double
eccentricity_with_symmetries_multithread(lbio_size_t n, std::string alphabet, lbio_size_t threads_) {
  double sum = 0;
  using PairVectorStrMult = std::vector<std::pair<std::string, lbio_size_t>>;
  using IterType = std::vector<std::pair<std::string, lbio_size_t>>::iterator;
  PairVectorStrMult tmp = permutation_invariant_strings_with_multiplicity(n, alphabet);
  PairVectorStrMult invariant_strings = remove_symmetry_invariant(tmp.begin(), tmp.end(), alphabet);
  std::vector<std::future<double>> v_futures;
  lbio_size_t M = invariant_strings.size();
  lbio_size_t m = static_cast<lbio_size_t>(std::ceil((double)M / (double)threads_));
  for (lbio_size_t i = 0; i < threads_-1; ++i) {
    auto it_begin = invariant_strings.begin() + m*i;
    auto it_end = invariant_strings.begin() + m*(i+1);
    v_futures.push_back(
      std::future<double>(std::async(
			    eccentricity_with_symmetries_iterator<IterType>, it_begin, it_end, n, alphabet)));
  }
  sum += eccentricity_with_symmetries_iterator(invariant_strings.begin() + m*(threads_-1),
					       invariant_strings.end(), n, alphabet);
  for (auto it = v_futures.begin(); it != v_futures.end(); ++it) {
    sum += it->get();
  }
  return sum / std::pow(alphabet.size(), n);
}

double
eccentricity_with_symmetries(lbio_size_t n, std::string alphabet, lbio_size_t threads_) {
  if (threads_ > 1) {
    return eccentricity_with_symmetries_multithread(n, alphabet, threads_);
  }
  auto tmp = permutation_invariant_strings_with_multiplicity(n, alphabet);
  auto invariant_strings = remove_symmetry_invariant(tmp.begin(), tmp.end(), alphabet);
  double sum = eccentricity_with_symmetries_iterator(invariant_strings.begin(), invariant_strings.end(),
						     n, alphabet);
  return sum / std::pow(alphabet.size(), n);
}



////////////////////////////////////////////////////////////////////////////////

// Useful alias used throughout the code
using ExactAlg    = EditDistanceWF<lbio_size_t, std::string>;
using BandApprAlg = EditDistanceBandApproxLinSpace<lbio_size_t, std::string>;

void
compare_edit_distance_algorithms(lbio_size_t n, lbio_size_t m, lbio_size_t k,
				 std::string alphabet, std::ostream& os) {
  lbio_size_t T_max = n / 2;
  lbio_size_t T_min = 1;

  ExactAlg exactAlg { n, m, {1,1,1} };

  GeometricProgression<lbio_size_t> geom(2, T_min);
  std::vector<lbio_size_t> Ts = geom.valuesLeq(T_max);
  Ts.push_back(0);
  

  
  lbio_size_t** dpMatrix = allocMatrix<lbio_size_t>(n+1, m+1); // !!!

  
  std::vector< std::shared_ptr<AlgorithmComparisonResult> > results;
  std::string s1(n, 'N');
  std::string s2(m, 'N');
  for (size_t l = 0; l < k; ++l) {
    std::shared_ptr<AlgorithmComparisonResult> res =
      std::make_shared<AlgorithmComparisonResult>();
    generateIIDString(s1, alphabet);
    generateIIDString(s2, alphabet);

    EditDistanceInfo tmp {};
    
    exactAlg.calculate(s1, s2);
    res->addExact(exactAlg.backtrack());

    // Approximation for all values of T
    for (auto T : Ts) {
      editDistanceBandwiseApproxMat(s1, s2, T, dpMatrix);
      closest_to_diagonal_backtrack(s1.size(), s2.size(), dpMatrix, tmp);
      res->addBandApprox(tmp, T);
    }
    results.push_back(res);
  }

  os << n << "\t";
  for (auto T : Ts) {
    os << T << "\t";
  }
  os << std::endl;
  for (auto pRes : results) {
    os << pRes->getExact() << "\t";
    for (auto T : Ts) {
      os << pRes->getBandApproxWithT(T) << "\t";
    }
    os << std::endl;
  }
  
  freeMatrix<lbio_size_t>(n+1, m+1, dpMatrix);
}


void
closest_to_diagonal_backtrack(size_t n, size_t m, size_t** dpMatrix,
			   EditDistanceInfo& info) {
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

      
lbio_size_t
optimal_bandwidth_exact(lbio_size_t n, std::string alphabet, double precision, lbio_size_t Tmin) {
  lbio_size_t T = Tmin;
  lbio_size_t T_2 = static_cast<lbio_size_t>(std::floor(n / 2));
  BandApprAlg exactAlg {n, n, T_2 ,{1,1,1}};
  IidPairGenerator gen(n, n, alphabet);
  lbio_size_t k_min = 5;
  while (T < n / 2) {
    double avg = 0;

    BandApprAlg apprAlg { n, n, T, {1,1,1} };
    for (lbio_size_t k = 0; k < k_min; ++k) {
      auto strings = gen();

      lbio_size_t exact = exactAlg.calculate(strings.first, strings.second);
      lbio_size_t approx = apprAlg.calculate(strings.first, strings.second);
      avg += (approx - exact) / static_cast<double>(exact);
    }
    if ( avg < k_min * precision) {
      return T;
    }
    T *= 2;
  }
  return T;
}

lbio_size_t
optimal_bandwidth(lbio_size_t n, std::string alphabet, double precision, lbio_size_t k, lbio_size_t Tmin) {
  lbio_size_t T_2 = Tmin;
  lbio_size_t T = 2 * T_2;
  
  BandApprAlg alg_T_2(n, n, T_2, {1,1,1});
  BandApprAlg alg_T(n, n, T, {1,1,1});
  IidPairGenerator gen(n,n,alphabet);

  while (T < static_cast<lbio_size_t>(std::floor(n/2.0))) {
    double avg = 0;
    for (lbio_size_t i = 0; i < k; ++i) {
      auto strings = gen();
      lbio_size_t ed_T = alg_T.calculate(strings.first, strings.second);
      lbio_size_t ed_T_2 = alg_T_2.calculate(strings.first, strings.second);
      avg += static_cast<double>(ed_T_2 - ed_T) / static_cast<double>(ed_T_2);
    }
    avg /= static_cast<double>(k);
    if (avg < precision) {
      return T_2;
    }
    T_2 = T;
    T *= 2;
    
  }
    
  
  return T_2;
}

     
} } } // namespaces
