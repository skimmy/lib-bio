#ifndef SIM_EDIT_H
#define SIM_EDIT_H

#include <memory>
#include <iostream>
#include <algorithm>

// Edit distance flgas
#define EDIT_DISTANCE_ESTIMATE_EXHAUSTIVE 0x1
#define EDIT_DISTANCE_ALGORITHM_QUADRATIC 0x2
#define EDIT_DISTANCE_INFO_PARTIAL        0x4
#define EDIT_DISTANCE_INFO_SCRIPT         0x8
#define EDIT_DISTANCE_SAMPLE_MATRIX       0x10 // 16

#define EDIT_DISTANCE_BOUNDED_ERROR       0x20 // 32
#define EDIT_DISTANCE_DIFF_BOUNDED_ERROR  0x40 // 64

// edit distance subtasks
#define EDIT_DISTANCE_SUBTASK_DEFUALT       0
#define EDIT_DISTANCE_SUBTASK_SCRIPT_DIST   8
#define EDIT_DISTANCE_SUBTASK_COMPARE_ALGS  32


/**
 * This is a structure to store information about how edit distance is
 * divided into substitution,
 */
class EditDistanceInfo
{
public:
  size_t n_sub;
  size_t n_del;
  size_t n_ins;  

  EditDistanceInfo() : n_sub(0), n_del(0), n_ins(0) {}
  EditDistanceInfo(const EditDistanceInfo& i) :
    n_sub(i.n_sub), n_del(i.n_del), n_ins(i.n_ins), edit_script(i.edit_script) {}

  std::string edit_script = "";

  size_t distance() { return n_sub + n_ins + n_del; }
  void reset() { n_sub = 0; n_del = 0; n_ins = 0; }
  bool operator==(const EditDistanceInfo& i) {
    return (n_sub == i.n_sub && n_del == i.n_del && n_ins == i.n_ins);
  }
  bool operator!=(const EditDistanceInfo& i) {
    return !(*this == i);
  }
  EditDistanceInfo& operator+=(const EditDistanceInfo& rhs) {
    n_sub += rhs.n_sub;
    n_del += rhs.n_del;
    n_ins += rhs.n_ins;
    return *this;
  }

  friend EditDistanceInfo operator+(EditDistanceInfo lhs, const EditDistanceInfo& rhs) {
    lhs += rhs;
    return lhs;
  }

 
  friend std::ostream& operator<<(std::ostream& out, const EditDistanceInfo& info) {
    out << info.n_sub << " " << info.n_del << " " << info.n_ins;
    return out;
  }  
};

//////////////////////////////////////////////////////////////////////
//
//             EDIT DISTANCE STRUCTS AND CLASSES
//
//////////////////////////////////////////////////////////////////////

// This is an attempt to make different edit distance algorithms fit a
// common interface. This is a work-in-progresso and is subject to
// many changes. The ultimate goal is to obtain an infrastructure that
// allows (almost) seamingless change of the algorithm without changes
// to the client code.

// Define macros for the indexes of costs in the cost vector. Used for
// mnemonics. This is not a good design and should be andanoned once a
// better way to pass costs will b eimplemented
#define iS_ 0
#define iD_ 1
#define iI_ 2

/** 
 * \brief This struct contains the minimal amount of data to support
 * dynamic programming algorithms, namely: the matrix and its sizes.
*/
template<typename T>
struct DynamicProgramming {
  T** dp_matrix;
  lbio_size_t n;
  lbio_size_t m;

  DynamicProgramming(lbio_size_t n_, lbio_size_t m_)
    : n(n_), m(m_)
  {
    dp_matrix = allocMatrix<T>(n+1,m+1);
  }

  ~DynamicProgramming() { freeMatrix<T>(n, m, dp_matrix); }
};

/**
 * \brief This class represents the standard Wagner and Fischer edit
 * distance dynamic programming algorithm.
 */
template<typename CostType, typename IndexedType>
class EditDistanceWF {
public:
    typedef std::vector<CostType> CostVector;

private:
  DynamicProgramming<CostType> dp_struct;
  // costs are in vector [W_S, W_D, W_I]
  // (i.e., [0] -> Sub, [1] -> Del, [2] -> Ins).
  // We should consider a better design for costs possible alternatives are
  // - a custom struct with fields for the csos
  // - a caller that is invoked when costs are needed (should be constexpr)
  CostVector costs_vector;
  
public:
  
  EditDistanceWF(lbio_size_t n, lbio_size_t m)
    : dp_struct {n, m}, costs_vector {1, 1, 1}
  {
    init();
  }

  EditDistanceWF(lbio_size_t n, lbio_size_t m, const CostVector& costs)
    : dp_struct{n, m}, costs_vector {costs}
  {
    init();
  }

  
  void init() {
    dp_struct.dp_matrix[0][0] = 0;
    for (lbio_size_t i = 1; i <= dp_struct.n; ++i) {
      dp_struct.dp_matrix[i][0]
	= dp_struct.dp_matrix[i-1][0] + costs_vector[iD_];
    }
    for (lbio_size_t j = 1; j <= dp_struct.m; ++j) {
      dp_struct.dp_matrix[0][j]
	= dp_struct.dp_matrix[0][j-1] + costs_vector[iI_];
    }
  }

  CostType calculate(const IndexedType& s1, const IndexedType& s2) {
    lbio_size_t n = s1.size();
    lbio_size_t m = s2.size();
    for (lbio_size_t i = 1; i <= n; ++i) {
      for(lbio_size_t j = 1; j <= m; ++j) {
	CostType delta = ( s1[i-1] == s2[j-1] ) ? 0 : costs_vector[iS_];
	CostType A_ = dp_struct.dp_matrix[i-1][j-1] + delta; 
	CostType B_ = dp_struct.dp_matrix[i-1][j] + costs_vector[iD_];
	CostType C_ = dp_struct.dp_matrix[i][j-1] + costs_vector[iI_];
	dp_struct.dp_matrix[i][j] = std::min(A_,std::min(B_, C_ ));
      }
    }
    return dp_struct.dp_matrix[n][m];
  }

  void print_dp_matrix() {
    printMatrix<CostType>(dp_struct.n+1, dp_struct.m+1, dp_struct.dp_matrix);
  }
}; // end of EditDistanceWF

template<typename CostType, typename IndexedType>
class EditDistanceBandApproxLinSpace {
public:
  typedef std::vector<CostType> CostVector;
  
private:
  DynamicProgramming<CostType> dp_struct;
  CostVector costs;
  lbio_size_t bandwidth;
  lbio_size_t n;
  lbio_size_t m;
  const CostType Inf;

public:
  EditDistanceBandApproxLinSpace(lbio_size_t n_, lbio_size_t m_, lbio_size_t T)
    : dp_struct {2, std::max(n_,m_)}, costs {1,1,1}, bandwidth {T}, n {n_}, m {m_},
      Inf {(*std::max_element(costs.begin(), costs.end())) * 2*(n_+m_+1)}
  {
  }

  void init() {    
    dp_struct.dp_matrix[0][0] = 0;
    for (lbio_size_t j = 1; j <= bandwidth; ++j) {
      dp_struct.dp_matrix[0][j] = dp_struct.dp_matrix[0][j-1] + costs[iI_];
    }
    /*for (lbio_size_t j = bandwidth+1; j <= m; ++j) {
      dp_struct.dp_matrix[0][j] = Inf;
      }*/
  }

  CostType calculate(const IndexedType& s1, const IndexedType& s2) {
    n = s1.size();
    m = s2.size();
    
    
    // This is a design decision but since at each calculation the DP
    // matrix (which are two vectors in this case) needs
    // initialization (which is not needed for the entire matrix
    // approach), the calculate method performs init at a cost of
    // wasting some time. For this reason caller may skip the init and
    // leave the calculate method to take care of initialization
    init();
    for (lbio_size_t i = 1; i <= n; ++i) {
      // 'Standard' initialization of the first column 
      dp_struct.dp_matrix[1][0] = (i <= bandwidth) ?
	dp_struct.dp_matrix[0][0] + costs[iD_] : Inf;

      // This takes care of 'diagonals' initialization (upper and
      // lower respectively)

      // TODO: Check the negative case
      if (i <= m - bandwidth) {
	dp_struct.dp_matrix[0][i+bandwidth] = Inf;
      }

      if (i >= bandwidth+1) {
	// Second index should always be >= 0 because of the if
	dp_struct.dp_matrix[1][i-(bandwidth+1)] = Inf;
	
      }
	
	
      // the fact the we used unsigned for indexes here becomes a pain
      // but I'd rather explicitly cast them back to int when
      // performed substraction than changing semantic. The j_max
      // calculation does not suffer of this problem because result
      // can not be negative (there are no differences involved)
      lbio_size_t j_min = std::max<int>(1, static_cast<int>(i - bandwidth));
      lbio_size_t j_max = std::min<int>(m, i + bandwidth);
      for (lbio_size_t j = j_min; j <= j_max; ++j) {
	CostType delta { (s1[i-1] == s2[j-1]) ? 0 : costs[iS_] };
	dp_struct.dp_matrix[1][j]
	  = std::min<CostType>(delta + dp_struct.dp_matrix[0][j-1],
			       std::min(dp_struct.dp_matrix[0][j] + costs[iD_],
					dp_struct.dp_matrix[1][j-1] + costs[iI_]));
      }
      
      // swap the two vectors
      std::swap<CostType*>(dp_struct.dp_matrix[0], dp_struct.dp_matrix[1]);
    }
    
    // TODO: Calculation
    return dp_struct.dp_matrix[0][m];
  }

  void print_dp_matrix() const {
    printMatrix<CostType>(2, m + 1, dp_struct.dp_matrix);
  }

  // This is here mostly for debugging purpose so it is not optimized
  // and may be removed in the future
  void resetMatrix() {
    for (lbio_size_t j = 0; j <= m; ++j) {
      dp_struct.dp_matrix[0][j] = 0;
      dp_struct.dp_matrix[1][j] = 0;
    }    
  }
			       
}; // end of  EditDistanceBandApproxLinSpace 

/**
 * \brief Extracts the substitutions array from an array of EditDistanceInfo
 */
std::unique_ptr<double[]> extractSubstitutionArray(const EditDistanceInfo* v, size_t k);

/**
 * \brief Extracts the deletions array from an array of EditDistanceInfo
 */
std::unique_ptr<double[]> extractDeletionArray(const EditDistanceInfo* v, size_t k);

/**
 * \brief Extracts the insertions array from an array of EditDistanceInfo
 */
std::unique_ptr<double[]> extractInsertionArray(const EditDistanceInfo* v, size_t k);


/**
 * \brief computes the edit distance between strings s1 and s2
 */
size_t editDistance(const std::string& s1, const std::string& s2);

size_t editDistanceEncoded(uint64_t s1, size_t n1, uint64_t s2, size_t n2, size_t** dpMatrix);

size_t editDistanceLinSpace(const std::string& s1, const std::string& s2, size_t* v0, size_t* v1);

EditDistanceInfo editDistanceLinSpaceInfo(const std::string& s1, const std::string& s2, EditDistanceInfo* v0, EditDistanceInfo* v1, EditDistanceInfo** sampleMat = NULL);

/**
 * \brief Computes the dynamic programming matrix dpMatrix for edit distance
 * between s1 and s2
 */
void editDistanceMat(const std::string& s1, const std::string& s2, size_t** dpMatrix);

void editDistanceWithInfo(const std::string& s1, const std::string& s2, EditDistanceInfo& info);

/**
 * \breif Reconstruct the edit script from the dynamic programming matrix.
 * The script and other informations are stored in the passed EditDistanceInfo 
 * structure.
 *
 * \param dpMatrix the dynamic programming matrix
 * \param n number of rows of the matrix minus one 
 * \param m number of columns of the matrix minus one
 * \param info the structure that will be filled with the information
 * 
 */
void
editDistanceBacktrack(size_t** dpMatrix, const std::string& s1, const std::string& s2,
		      EditDistanceInfo& info);

void
closestToDiagonalBacktrack(size_t n, size_t m, size_t** dpMatrix, EditDistanceInfo& info);

size_t
editDistanceBandwiseApprox(const std::string& s1, const std::string& s2, size_t T);


double
testExhaustiveEditDistanceEncoded(size_t n, double* freq);

void
computeAverageDPMatrix(double** dpMatrix, size_t n, size_t m);

SampleEstimates
editDistanceErrorBoundedEstimates(size_t n, double precision, double z_delta, size_t k_min = 16);

SampleEstimates
editDistanceRelativeErrorEstimates(size_t n, double e_model, double precision, double z_delta);

/**
 * Computes a matrix containing in (i,j) the number of times a
 * 'closest diagonal' path traverses the cell (i,j)
 *
 * @param n number of rows of the matrix
 * @param m number of columns of the matrix
 * @param k number of samples to be used
 * @param distMatrix pointer to a matrix that will contain the frequencies (i.e., output)
 * @param scripts a pointer to a standard vector that will contain the scripts. 
 *  If set to <code>nullptr</code> no script will be stored
 */
void
scriptDistributionMatrix(size_t n, size_t m, size_t k, size_t** distMatrix,
			 std::vector<std::string>* scripts = nullptr);


void
compareEditDistanceAlgorithms(size_t n, size_t m, size_t k, std::ostream& os = std::cout);

//////////////////////////////////////////////////////////////////////
//            EDIT DISTANCE ESTIMATION AND SAMPLING
//////////////////////////////////////////////////////////////////////

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
void
editDistanceEstimations(size_t n_min, size_t n_max, size_t n_step, size_t k_max);

/*
 * \brief Computes k_the edit distance for k_samples pairs of random strings
 * with length n. The distances are stored in an array whose (smart) pointer 
 * is returned.
 *
 * \param n the length of the strings
 * \param k_samples the number of pairs to compute
 * \return a (smart) pointer to an array containing all the distances
 */
std::unique_ptr<size_t[]>
editDistSamples(size_t n, size_t k_samples);

std::unique_ptr<EditDistanceInfo[]>
editDistSamplesInfo(size_t n, size_t k_samples);

std::unique_ptr<EditDistanceInfo[]>
editDistSamplesInfoLinSpace(size_t n, size_t k_samples, EditDistanceInfo** sampleMat = NULL);

/**
 * \brief This class is a generator of edit distances, it can be used
 * to sample the edit distance distribution for gven n and m. 

 * The class is not to be used when also partition of the operation is
 * needed. 

 * The class is templated by an algorithm that actually computes the
 * distance. This class must be an instantiation of algorithms class
 * defined above with cost type `lbio_size_t` and indexed type
 * `std::string` moreover.
 * - the class must have a constructor accepting dimensions n and m of type
 *   `lbio_size_t` (or type implicitly converting from `lbio_size_t`)
 * - it must contain a method `calculate` accepting two indexed input
 */

template<typename EDAlg_>
class EditDistanceSample {
private:
  std::string A_;
  std::string B_;
  
public:
  EditDistanceSample(lbio_size_t n, lbio_size_t m)
    : A_(n, 'N'), B_(m, 'N') {  }
  
  lbio_size_t operator()(EDAlg_& algorithm) {
    generateIIDString(A_);
    generateIIDString(B_);
    return algorithm.calculate(A_, B_);
  }

  std::pair<std::string, std::string> latest_strings() {
    return std::pair<std::string, std::string>(A_, B_);
  }
  
};

/**
 * \brief Convience struct to manage vecotrs, this is a bad design and
 * should be changed.
 */
struct EditDistanceSimOutput {
  double* distPDF = NULL;
  
  ~EditDistanceSimOutput() {
    delete[] this->distPDF;
    this->distPDF = NULL;
  }
};

#endif
