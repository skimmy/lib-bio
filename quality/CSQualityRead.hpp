#ifndef CS_QUALITY_READ_H
#define CS_QUALITY_READ_H

#include "ReadQuality.hpp"

#include <string>
using namespace std;

/**
 * \deprecated{see Quality and PhredQuality}
 * \brief This class represents a quality sequence expressed with Phred values
 * (\b deprecated see Quality and PhredQuality)
 *
 * \sa Quality
 * \sa PhredQuality
 */
class CSQualityRead : public ReadQuality {
public:
  // ---------------------------------------------------------
  //                       CONSTRUCTORS
  // ---------------------------------------------------------
  /**
   * \brief Constructs a CSQualityRead for a \c n long sequence
   * from the given quality values tring
   *
   * Note that the passed quality values string must contain at
   * least \c n qualitiy values otherwise unpredictable behavior
   * may arise.
   *
   * \param n The length of the sequence
   * \param qualities The quality values string
   *
   * \sa CSQualityRead(const CSQualityRead& other)
   */
  CSQualityRead(size_t n, const string& qualities);
  /**
   * \brief Creates a copy of a given CSQualityRead object
   *
   * \param other The CSQualityRead to be copied
   *
   * \sa CSQualityRead(size_t n, const string& qualities)
   */
  CSQualityRead(const CSQualityRead& other);
private:
  // ---------------------------------------------------------
  //                  UTILITY PRIVATE METHODS
  // ---------------------------------------------------------
  void computeValues(const string& quals);
  void computeProbabilities();
};

#endif
