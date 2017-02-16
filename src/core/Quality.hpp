#ifndef QUALITY_CLASS_H
#define QUALITY_CLASS_H

#include <string>

/**
 * \brief This is a virtual class (interface) defining the
 * common operations that can be performed with a sequence of
 * qualities (regardless the specific representation).
 */
class Quality {
public:
  virtual ~Quality() {}
  /**
   * \brief Returns a vector containing the error probabilities.
   *
   * This methods can be used with zero, one or two parameters.
   * When called without any parameter it returns the entire
   * stored vector converted into a \em probabilistic encoding,
   * when used with one paramater programmer can decide from which
   * position the vector should start (spanning until the end
   * of the sequence), when called with two parameters it returns
   * a vector containing the specified number of elements filled
   * with the corresponding probabilities. If the size of returned
   * vector is greater than the total number of position available
   * the remaining elements <b>are not initialized</b>.
   *
   * If n is the number of elements stored in the Quality object
   * the three following calls are equivalent (and all return
   * the entire sequence converted into a probabilities vector).
   * \code
   * getProbabilities();
   * getProbabilities(0);
   * getProbabilities(0,n);
   * \endcode
   *
   * \param begin The starting position from which convert the
   * stored values
   * \param length The length of the returned vector (which is
   * the length if internal sequence if it is no specified)
   * \return The conversion into the probabilities representation
   * of the internal sequence.
   */
  virtual double* getProbabilities(size_t begin = 0, size_t length = 0) const = 0;

  /**
   * \brief Returns a vector containing the quality values.
   *
   * Thid method can be used with zero, one or two parameters
   * When called without any parameter it returns the entire
   * stored vector converted into a \em quality encoding when
   * used with one paramater programmer can decide from which
   * position the vector should start (spanning until the end
   * of the sequence), when called with two parameters it returns
   * a vector containing the specified number of elements filled
   * with the corresponding qualities . If the size of returned
   * vector is greater than the total number of position available
   * the remaining elements <b>are not initialized</b>.
   *
   * If n is the number of elements stored in the Quality object
   * the three following calls are equivalent (and all return
   * the entire sequence converted into a probabilities vector).
   * \code
   * getQualities();
   * getQualities(0);
   * getQualities(0,n);
   * \endcode
   *
   * \param begin The starting position from which convert the
   * stored values
   * \param length The length of the returned vector (which is
   * the length if internal sequence if it is no specified)
   * \return The conversion into the qualities representation
   * of the internal sequence.
   */
  virtual int* getQualities(size_t begin = 0, size_t length = 0) const = 0;

  /**
   * \brief Returns the overall probability for the entire sequence.
   *
   * The overall probability depends, of course, by the model
   * used for the error probability, we will assume (until
   * otherwise stated) that the model is an <em>indipendend
   * probabilities model</em> and, therefore, the returned value
   * is the product of all the values (interpreted as probabilities)
   * stored in the internal sequence.
   *
   * \return The overall probability
   */
  virtual double getOverallProbability() = 0;
  /**
   * \brief Returns the length of the internal sequence.
   *
   * \return The length of the internal sequence
   */
  virtual size_t length() const = 0;
};

#endif
