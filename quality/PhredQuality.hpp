#ifndef PHRED_QUALITY_H
#define PHRED_QUALITY_H

#include "Quality.hpp"

#include <cstring>
#include <string>
using namespace std;

/**
 * \brief This class represents the qualty of a sequence using
 * the \em Phred coding.
 *
 * This class should be used when natively the quality is threated
 * a a \em Phred encoded sequence because the internal repreentation
 * is given in phred integer values. Conversion to a probabilistic
 * representation is therefore a little inefficient since it
 * requires computing probabilistic valeus starting from the
 * integer representation.
 * The conversion is performed using the phred to probability
 * equation 
 * \f[ p = 10^{-\frac{q}{10}} \f]
 * and therefore the probability should be interpreted as <b>error
 * probability</b>
 *
 * \sa Quality
 * \sa ProbabilisticQuality
 */
class PhredQuality : public Quality {
private:
  int* qualVector;
  size_t n;

public:
  // ---------------------------------------------------------
  //                 CONSTRUCTOS AND DESTRUCTOR
  // ---------------------------------------------------------

  /**
   * \brief This constructor creates a PhredQuality object starting
   * from an integer vector.
   *
   * \param quals The vector of quality values
   * \param n The size of the quality values vector
   *
   * \sa PhredQuality(const string& quals, size_t n)
   * \sa PhredQuality(const PhredQuality& other)
   */
  PhredQuality(const int* quals, size_t n);
  /**
   * \brief This constructor creates a PhredQuality object from a
   * string containing the quality values. 
   *
   * The input string must contain a number of integer as indicated in the
   * second parameter, numbers must be separated using space
   * characters (i.e. blank, tab, newlen , ... ).
   *
   * \param quals The string containing quality values
   * \param n The number of values containing in the input string
   * 
   * \sa PhredQuality(const int* quals, size_t n)
   * \sa PhredQuality(const PhredQuality& other)
   */
  PhredQuality(const string& quals, size_t n);
  /**
   * \brief The copy constructor creates a copy of the passed PhredQuality
   * instance.
   *
   * \param other The PhredQuality to be copied
   *
   * \sa PhredQuality(const int* quals, size_t n)
   * \sa PhredQuality(const string& quals, size_t n);
   */
  PhredQuality(const PhredQuality& other);
  ~PhredQuality();

  // ---------------------------------------------------------
  //                 'QUALITY' OVERRIDE METHODS
  // ---------------------------------------------------------
  double* getProbabilities(size_t begin = 0, size_t length = 0) const;
  int* getQualities(size_t begin = 0, size_t length = 0) const;
  double getOverallProbability() const;
  size_t length() const;

  // ---------------------------------------------------------
  //                       STATIC METHODS
  // ---------------------------------------------------------

  /**
   * \brief Converts a vector represented as a probability value into
   * a vector representing phred encoded qualities.
   *
   * Note that the returned vector is dynamically allocated and should,
   * therefore, be deleted by the programmer in order to prevent
   * memory leaks.
   *
   * \param probs The vector containing the probabilities
   * \param n The size of the input (and output) vector
   * \return A vector of integers containing the conversion
   * into phred values of the input vector
   */
  static int* toPhred(const double* probs, size_t n);
  /**
   * \brief Converts a single quality value expressed as a probability
   * into an integer value representing the phred encoded quality.
   *
   * \param p The probability to be converted
   * \return The phred quality value if the input
   */
  static int toPhred(double p);

private:
  // ---------------------------------------------------------
  //                      UTILITY METHODS
  // ---------------------------------------------------------
  void init(const int* v, size_t n);
  void parseQualityString(const string& s);
};

#endif
