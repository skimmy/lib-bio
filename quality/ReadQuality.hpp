#ifndef READ_QUALITY_H
#define READ_QUALITY_H

#include <cstdint>
#include <stdlib.h>

extern double probLookup[256];
void initLookupTables();

/**
 * \deprecated{Use Quality}
 * \brief This class Represents the quality of read (\b Deprecated use Quality instead)
 *
 * \sa Quality
 * \sa PhredQuality
 * ProbabilisticQuality
 */
class ReadQuality {
protected:
  /**
   * \brief The quality values vector
   */
  uint8_t* qualities;
  /**
   * \brief The probability values vector
   */
  double* probabilities;
  /**
   * \brief The length of the quality sequence
   */
  size_t length;
public:
  // ---------------------------------------------------------
  //                 CONSTRUCTORS AND DESTRUCTOR
  // ---------------------------------------------------------
  /**
   * \brief Createsa an empty ReadQuality object with given size
   *
   * \param n The size of the quality sequence
   *
   * \sa ReadQuality(const ReadQuality& other)
   */
  ReadQuality(size_t n);
  /**
   * \brief Creates a copy of the given ReadQuality
   *
   * \param other The ReadQuality to be copied
   *
   * \sa ReadQuality(size_t n)
   */
  ReadQuality(const ReadQuality& other);
  ~ReadQuality();
  
  // ---------------------------------------------------------
  //                    GET AND SET METHODS
  // ---------------------------------------------------------
  /**
   * \brief Returns the length of the quality sequence
   *
   * \return The length of the quality sequence
   *
   * \sa getQualities()
   * \sa getProbabilities()
   */
  size_t getLength() const;
  /**
   * \brief Returns the quality values sequence
   *
   * \return The quality values sequence
   * 
   * \sa getProbabilities()
   * \sa getLength()
   */
  uint8_t* getQualities();
  /**
   * \brief Returns the probabilities sequence
   *
   * \return The probabilities sequence
   *
   * \sa getQualities()
   * \sa getLength()
   */
  double* getProbabilities();
  /**
   * \brief Returns the overall probability of the quality sequence
   *
   * \return The overall probability
   *
   * \sa getProbabilities()
   * \sa getQualities()
   */
  double getOverallProbability() const;
protected:
  // ---------------------------------------------------------
  //                      UTILITY METHODS
  // ---------------------------------------------------------
  /**
   * \brief Initializes the internal structure of the quality
   * sequence
   *
   * \param n The length of the internal structure
   */
  void init(size_t n);
};

#endif
