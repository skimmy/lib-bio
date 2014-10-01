#ifndef MATCH_SIMILARITY_H
#define MATCH_SIMILARITY_H

#include "../adt.h"

#include <string>
using namespace std;

/**
 * \brief This class represents a <em>match similarity</em>
 * dynamic programming algorithm
 *
 * Match similiraty is the numbero fo mismatches between two
 * sequences. This class implements the dynamic programming
 * algorithm to compute this \e metric.
 * The class inherites from DynamicProgramming class in order
 * to provide the generic <em> dynamic programming</em> interface.
 *
 * \sa DynamicProgramming
 * \sa SmithWatermanDP
 */
class MatchSimilarity : public DynamicProgramming {
private:
  int** matrix;
  string* seq1;
  string* seq2;
public:
  // ---------------------------------------------------------
  //                  CONSTRUCTORS DESTRUCTOR
  // ---------------------------------------------------------
  /**
   * \brief Creates the MatchSimilarity object from the two
   * sequences given.
   *
   * For efficiency reasons the two sequences are not locally
   * copied and therefore the passed references should reamin
   * valid during the entire execution of the algorithm (i.e
   * the object scope).
   *
   * \param s1 The first sequence
   * \param s2 The second sequence
   *
   */
  MatchSimilarity(string* s1, string* s2);
  ~MatchSimilarity();

  // ---------------------------------------------------------
  //            override from 'DynamicProgramming'
  // ---------------------------------------------------------
  void initMatrix();
  void computeEntry(size_t i, size_t j);
  void printMatrix();
};

#endif
