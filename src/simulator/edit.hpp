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

#define EDIT_DISTANCE_BOUNDED_ERROR       0x20 // 32
#define EDIT_DISTANCE_DIFF_BOUNDED_ERROR  0x40 // 64

#define EDIT_DISTANCE_BANDWIDTH_ESTIMATE  0x80 // 128
// Approx ~g(n) --> 64+128 = 192 --> -f 192

// edit distance subtasks
#define EDIT_DISTANCE_SUBTASK_DEFUALT       0
#define EDIT_DISTANCE_SUBTASK_SCRIPT_DIST   8
#define EDIT_DISTANCE_SUBTASK_COMPARE_ALGS  32

// !!! Forward declaration to be removed after namespace enclosing is done !!!!
class EditDistanceInfo;


// These namespace will eventually contain all functions and classes
namespace lbio { namespace sim { namespace edit {


/**
   \brief Performs backtrack and fills the \c EditDistanceInfo
   passed. When ties are possible the backtrack goes through the
   diagonal closest to the main one.

   \param n         Number of rows of matrix minus one
   \param m         Number of columns of matrix minus one
   \param dpMatrix  The dynamic programming matrix
   \param info      The object used to put output
 */
void
closest_to_diagonal_backtrack(size_t n, size_t m, size_t** dpMatrix,
			      EditDistanceInfo& info);



//////////////////////////////////////////////////////////////////////
//           COMPARISON OF ALGORITHMS FOR EDIT DISTANCE
//////////////////////////////////////////////////////////////////////

void
compareEditDistanceAlgorithms(size_t n, size_t m, size_t k, std::ostream& os = std::cout);


/**
   \brief Computes the first value T such that the difference between 
 */
lbio_size_t
optimal_bandwidth(lbio_size_t n, double precision, lbio_size_t Tmin = 1);


} } } //namespaces

/**
 * \brief This is a class to store information about how edit distance
 * is divided into substitution, deletions and insertions.
 *
 * The members are kept public since they are often accessed and this
 * used to be a struct.
 */
class EditDistanceInfo
{
public:
  lbio_size_t n_sub;
  lbio_size_t n_del;
  lbio_size_t n_ins;
  std::string edit_script;

  EditDistanceInfo();
  EditDistanceInfo(lbio_size_t s, lbio_size_t d, lbio_size_t i);
  EditDistanceInfo(lbio_size_t c);

  size_t
  distance() const;
  
  void
  reset();

  bool
  operator<(const EditDistanceInfo& i) const;
  
  bool
  operator==(const EditDistanceInfo& i) const;

  bool
  operator!=(const EditDistanceInfo& i) const;
  
  EditDistanceInfo&
  operator+=(const EditDistanceInfo& rhs);

  friend EditDistanceInfo
  operator+(EditDistanceInfo lhs, const EditDistanceInfo& rhs);

  EditDistanceInfo&
  operator*=(lbio_size_t scalar);

  friend EditDistanceInfo
  operator*(EditDistanceInfo lhs, lbio_size_t rhs);
 
  friend std::ostream&
  operator<<(std::ostream& out, const EditDistanceInfo& info);
};

const EditDistanceInfo InfoUnitSub {1,0,0};
const EditDistanceInfo InfoUnitDel {0,1,0};
const EditDistanceInfo InfoUnitIns {0,0,1};

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
 *
 * The members are public so that matrix and its sizes can be
 * acccessed and using directly. It is however reccomended to use the
 * methods and operator to guarantee maximum portability of the
 * software (e.g., in case of internal changes to the struct)
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

  /**
   * \brief This operator allows accessing (read and write) elements
   * of the dynamic programming matrix.
   */
  T& operator()(lbio_size_t i, lbio_size_t j) {
    return dp_matrix[i][j];
  }

  void swap_rows(lbio_size_t i1, lbio_size_t i2) {
    std::swap<T*>(dp_matrix[0], dp_matrix[1]); 
  }

  void print_matrix() {
    printMatrix<T>(2, m + 1, dp_matrix);
  }
}; // DynamicProgramming D



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

  EditDistanceWF(lbio_size_t n, lbio_size_t m, const CostVector& costs)
    : dp_struct{n, m}, costs_vector {costs}
  {
    init();
  }

  
  void
  init() {
    dp_struct(0, 0) = 0;
    for (lbio_size_t i = 1; i <= dp_struct.n; ++i) {
      dp_struct(i, 0)
	= dp_struct(i-1, 0) + costs_vector[iD_];
    }
    for (lbio_size_t j = 1; j <= dp_struct.m; ++j) {
      dp_struct(0, j)
	= dp_struct(0, j-1) + costs_vector[iI_];
    }
  }

  CostType
  calculate(const IndexedType& s1, const IndexedType& s2) {
    lbio_size_t n = s1.size();
    lbio_size_t m = s2.size();
    for (lbio_size_t i = 1; i <= n; ++i) {
      for(lbio_size_t j = 1; j <= m; ++j) {
	CostType delta = ( s1[i-1] == s2[j-1] ) ? 0 : costs_vector[iS_];
	CostType A_ = dp_struct(i-1, j-1) + delta; 
	CostType B_ = dp_struct(i-1, j)   + costs_vector[iD_];
	CostType C_ = dp_struct(i, j-1)   + costs_vector[iI_];
	dp_struct(i, j) = std::min(A_,std::min(B_, C_ ));
      }
    }
    return dp_struct(n, m);
  }

  EditDistanceInfo
  backtrack() {
    EditDistanceInfo info;
    // following uses an 'old' functions that accepts size_t** as DP matrix
    // it should be changed to accept a dp_struct
    lbio::sim::edit::closest_to_diagonal_backtrack(dp_struct.n, dp_struct.m,
						   dp_struct.dp_matrix, info);
    return info;
  }

  void
  print_dp_matrix() {
    dp_struct.print_matrix();
  }
}; // EditDistanceWF

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
  EditDistanceBandApproxLinSpace(lbio_size_t n_, lbio_size_t m_,
				 lbio_size_t T, const CostVector& costV)
    : dp_struct {2, std::max(n_,m_)}, costs {costV}, bandwidth {T},
      n {n_}, m {m_},
      Inf {(*std::max_element(costs.begin(), costs.end())) * 2*(n_+m_+1)}
  { }

  void init() {    
    dp_struct(0, 0) = 0;
    for (lbio_size_t j = 1; j <= bandwidth; ++j) {
      dp_struct(0, j) = dp_struct(0, j-1) + costs[iI_];
    }
  }

  CostType calculate(const IndexedType& s1, const IndexedType& s2) {
    n = s1.size();
    m = s2.size();    
    
    // initialization is done here rather than in the constructor
    // because needed at each calculation (i.e., vectors will contain
    // values from older calculations if any)
    init();
    
    for (lbio_size_t i = 1; i <= n; ++i) {
      // initialization of 'border' and 'diagonals'
      dp_struct(1, 0) = (i <= bandwidth) ?
	dp_struct(0, 0) + costs[iD_] : Inf;
      if (i <= m - bandwidth) {
	dp_struct(0, i+bandwidth) = Inf;
      }
      if (i >= bandwidth+1) {
	dp_struct(1, i-(bandwidth+1) ) = Inf;	
      }
	
      // casts are necessary to cope with negative 
      lbio_size_t j_min = std::max<int>(1, static_cast<int>(i - bandwidth));
      lbio_size_t j_max = std::min<int>(m, i + bandwidth);
      for (lbio_size_t j = j_min; j <= j_max; ++j) {
	CostType delta { (s1[i-1] == s2[j-1]) ? 0 : costs[iS_] };
	dp_struct(1, j)
	  = std::min<CostType>(delta + dp_struct(0, j-1),
			       std::min(dp_struct(0, j) + costs[iD_],
					dp_struct(1, j-1) + costs[iI_]));
      }      
      // swap the two vectors
      dp_struct.swap_rows(0, 1);      
    }    
    return dp_struct(0, m);
  }

  
  void print_dp_matrix() const {
    dp_struct.print_matrix();
  }
			       
}; // EditDistanceBandApproxLinSpace



//////////////////////////////////////////////////////////////////////
//                      CONVERSION FUNCTIONS
//////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////
//              EDIT DISTANCE CALCULATION FUNCTION
//////////////////////////////////////////////////////////////////////

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


size_t
editDistanceBandwiseApprox(const std::string& s1, const std::string& s2, size_t T);


//////////////////////////////////////////////////////////////////////
//                          BACKTRACKING
//////////////////////////////////////////////////////////////////////

/**
   \brief Reconstruct the edit script from the dynamic programming matrix.
   The script and other informations are stored in the passed EditDistanceInfo 
   structure.
  
   \param dpMatrix the dynamic programming matrix
   \param n number of rows of the matrix minus one 
   \param m number of columns of the matrix minus one
   \param info the structure that will be filled with the information
   
 */
void
editDistanceBacktrack(size_t** dpMatrix, const std::string& s1,
		      const std::string& s2, EditDistanceInfo& info);

//////////////////////////////////////////////////////////////////////
//            EDIT DISTANCE ESTIMATION AND SAMPLING
//                 (TO BE MOVED TO EDIT_EST.HPP)
//////////////////////////////////////////////////////////////////////


double
testExhaustiveEditDistanceEncoded(size_t n, double* freq);

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
scriptDistributionMatrix(size_t n, size_t m, size_t k,
			 std::vector<std::string>& scripts);

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
  lbio::sim::generator::IidPairGenerator gen;
  
public:
  EditDistanceSample(lbio_size_t n, lbio_size_t m)
    : gen(n,m) {  }
  
  lbio_size_t operator()(EDAlg_& algorithm) {
    auto next = gen();
    return algorithm.calculate(next.first, next.second);
  }

  std::pair<std::string, std::string> latest_strings() {
    return gen.last_pair();
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

