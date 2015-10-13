#ifndef SMITH_WATERMAN_DP
#define SMITH_WATERMAN_DP

#include "../algorithms.h"

#include <string>

/**
 * \brief This class represents the implementation of a <em>Smith
 * Waterman</em> dynamic programming algorithm.
 *
 * The algorithm is implemented using an internal \e int matrix
 * and maintaining two pointers to the sequences to be compared.
 *
 * \sa DynamicProgramming
 * \sa MatchSimilarity
 */
class SmithWatermanDP 
  : public DynamicProgramming 
{
private: 
  int** matrix;
  const char* x;
  const char* y;
  // algorithm parameters
  int** sim;
  int gapPenalty;
  // backtrack related members
  BacktrackOperation** btMatrix;
  bool btEnabled;
public:
  // ---------------------------------------------------------
  //                CONSTRUCTORS AND DESTRUCTOR
  // ---------------------------------------------------------
  /**
   * \brief Constructs the SmithWatermanDP object from the two
   * input sequences.
   *
   * \param s1 The first sequence
   * \param s2 The second sequence
   * \param n1 The length of the first sequence
   * \param n2 The length of the second sequence
   */
  SmithWatermanDP(const char* s1, size_t n1, const char* s2, size_t n2);

  SmithWatermanDP(const std::string& s1, const std::string& s2);
  ~SmithWatermanDP();

  // ---------------------------------------------------------
  //            OVERRIDE FROM 'DYNAMICPROGRAMMING'
  // ---------------------------------------------------------
  void initMatrix();
  void computeEntry(size_t i, size_t j);
  void computeMatrix();
  void printMatrix();

  // ---------------------------------------------------------
  //          SW ALGORITHM PARAMETERS SETTING METHODS
  // ---------------------------------------------------------
  /**
   * \brief Sets the <em>similarity matrix</em> for the
   * Smith Waterman algorithm
   * 
   * The similarity matrix is used to assign a \e score to
   * the comparison between two smbols of the alphabet.
   * The matrix should therefore have proper dimensions (namely
   * \f$ k \times k \f$ where \f$ k \f$ is the size of the
   * alphabet from which symbols are drawn.
   * The matrix is neither copied neither deallocated, its
   * management is therefore left to the programmer.
   *
   * \sa setGapPenalty()
   *
   */
  void setSimilarityMatrix(int** s);
  /**
   * \brief Sets the <em>gap penalty</em> for the Smith
   * Waterman algorithm
   *
   * The gap penalty is used to give a score (usually negative)
   * to each of the gaps in the sequence.
   * 
   * \sa setSimilarityMatrix()
   */
  void setGapPenalty(int d);
  // ---------------------------------------------------------
  //                RESEULT RETRIEVING METHODS
  // ---------------------------------------------------------
  /**
   * \brief Returns the position.
   *
   * If two or more positions have same best score one of them is returned
   * (usually the first one in a row major scan of the matrix, but this is not
   * guaranteed and can change).
   *
   * \return The position of the best computed score
   *
   * \sa getScoreAt()
   */
  MatrixPoint2D getGlobalBest() const;
  
  /**
   * \brief Returns the score at the given position
   *
   * \param The MatrixPoint2D for the position we want to query the score.  
   *
   * At the current implementation there is no check for the indices of the given
   * position, if they are not on a valid range the behaviour is umpredexted.
   *
   * \return The score at the given position.
   *
   * \sa getGlobalBest()
   */
  int getScoreAt(const MatrixPoint2D& p);

  // ---------------------------------------------------------
  //                   BACKTRACKING METHODS  
  // ---------------------------------------------------------
  void enableBacktrack();
  void disableBacktrack();
  bool isBacktrackEnabled() const;
  void printBacktrackMatrix() const;
  
private:
  // ---------------------------------------------------------
  //                  PRIVATE UTILITY METHODS
  // ---------------------------------------------------------
  // matrix creation and destructon
  void createMatrix();
  void destroyMatrix();
  // algorithm parameters helper(s)
  int** createDefaultSimilarityMatrix(size_t s);
  // querying
  void indicesOfMaxElement(size_t& i, size_t& j) const;
  // backtracking
  void initBtMatrix();
  void deleteBtMatrix();
  void resetBtMatrix();
  BacktrackOperation getBtForPosition(size_t i, size_t j) const;
};

#endif
