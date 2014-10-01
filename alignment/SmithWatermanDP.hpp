#ifndef SMITH_WATERMAN_DP
#define SMITH_WATERMAN_DP

#include "../adt.h"

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
public:
  // ---------------------------------------------------------
  //                constructors and destructor
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
private:
  // ---------------------------------------------------------
  //                  PRIVATE UTILITY METHODS
  // ---------------------------------------------------------
  void createMatrix();
  void destroyMatrix();
  int** createDefaultSimilarityMatrix(size_t s);
};

#endif
