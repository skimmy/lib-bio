#ifndef PROBABILISTIC_QUALITY_H
#define PROBABILISTIC_QUALITY_H

#include "../quality.h"

#include <string>


/**
 * \brief This class represents the quality of a string represented
 * as a <em>probabilistic</em> measure.
 *
 * This class should be used to represent a quality string as
 * probabilistic qualities. Converting to \e phred qualities
 * should avoided as much as possible because it requires an
 * array instatiation and filling.
 * The conversion is perfomed using the probability to phred
 * equation
 * \f[ q = -10\log_{10}{p} \f]
 * and therefore the probability are actually <em>error
 * probbilities</em>
 *
 * \sa Quality
 * \sa PhredQuality 
 */
class ProbabilisticQuality : public Quality {
private:
  double* probVector;
  size_t n;
public:
  // ---------------------------------------------------------
  //                CONSTRUCTORS AND DESTRUCTOR
  // ---------------------------------------------------------
  /**
   * \brief Creates a ProbabilisticQuality object from a
   * vecotr of \c double
   *
   * \param prob The probabilities vector
   * \param n The size of the vector
   *
   * \sa ProbabilisticQuality(const ProbabilisticQuality& other)
   */
  ProbabilisticQuality(const double* prob, size_t n);
  /**
   * \brief Creates a copy of the given ProbabilisticQuality
   *
   * \param other The ProbabilisticQuality to copied
   *
   * \sa ProbabilisticQuality(const double* prob, size_t n)
   */
  ProbabilisticQuality(const ProbabilisticQuality& other);
  ~ProbabilisticQuality();

  // ---------------------------------------------------------
  //                'QUALITY' OVERRIDE METHODS
  // ---------------------------------------------------------
  double* getProbabilities(size_t begin = 0, size_t length = 0) const;
  int* getQualities(size_t begin = 0, size_t length = 0) const;
  double getOverallProbability() const;
  size_t length() const;

  // ---------------------------------------------------------
  //                      STATIC METHODS
  // ---------------------------------------------------------
  /**
   * \brief Converts a vector of \e phred qualities value into a vector
   * of probabilistic values.
   *
   * Note that the returned reference points to a dynamically
   * allocated array and should, therefore, freed as soon as
   * it is no more needed
   *
   * \param quals The qualitiy values to be converted
   * \param n The size of the input (and output) vector
   * \return A vector of \c double containing the conversione
   * into probabilistic qualities of the input phred qualities
   * 
   * \sa toProbabilistic(int qual)
   */
  static double* toProbabilistic(const int* quals, size_t n);
  /**
   * \brief Converts a single \e phred quality into a probabilistic
   * quality.
   *
   * \param qual The \e phred quality value to be converted 
   * \return The probability quality converted from input
   *
   * \sa toProbabilistic(const int* quals, size_t n)
   */
  static double toProbabilistic(int qual);


  // Static Factory methods
  static ProbabilisticQuality fromEncodedQuality(const std::string& v, qual::QualityEncodingType enc);

private:
  // ---------------------------------------------------------
  //                      UTILITY METHODS
  // ---------------------------------------------------------
  void init(const double* v, size_t n);
};

typedef ProbabilisticQuality ProbQuality;


#endif
