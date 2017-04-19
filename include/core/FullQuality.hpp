#ifndef FULL_QUALITY_H
#define FULL_QUALITY_H

#include <core/Quality.hpp>
#include <core/QualifiedSequence.hpp>
#include <core/Sequence.h>

#include <string>

/**
 * \brief This class represent a sequence of quality values
 *
 * The FullQuality class represent a sequence of quality values
 * expressed as error probabilities, for a given position (i.e
 * from 0 to length - 1) there are S probabilities where S is
 * the size of the alphabet C.
 * The class C passed as template parameter must implement all
 * the require static methods: 
 * \code
 * length()
 * getChaR(size_t)
 * getIndex(char)
 * \endcode
 *
 * \sa FullyQualifiedSequence
 * \sa DNAAlphabet
 * \sa ColorAlphabet
 */
template<class C>
class FullQuality {
private:
  size_t symbolCount;
  double** probVector;
  size_t length;
public:

  // ---------------------------------------------------------
  //              - CONSTRUCTORS AND DESTRUCTOR              -
  // ---------------------------------------------------------


  /**
   * \brief Construct an emtpy FullQuality sequence
   *
   * Even if the sequence itself is \em empty, the internal
   * state is initialized with in order to store a sequence
   * whose length is indicated by the parameter.
   *
   * \param n The length of the internal sequence 
   * 
   * \sa FullQuality(const Quality& qual, const Sequence& seq)
   * \sa FullQuality(const QualifiedSequence& qualSeq)
   */
  FullQuality(size_t n);
  /**
   * \brief construct the Fullquality from a Quality and a
   * Sequence object
   *
   * This constructor takes a Quality object and a Seqeunce as
   * parameter and creates the FullQuality from them. The way
   * probabilities are assigned is as follow:
   * - To the character `c` at the ith position of the sequence
   *   is given the corresponding ith probabilities \f$p\f$
   * - The \f$1 - p\f$ is then evenly partitioned among all but 
   *   `c` characters in the alphabet C
   * Note that is programmer responsability to ensure that the
   * alphabet used for the Sequence parameter is compatible with
   * the template parameter C.
   *
   * \param qual The Quality class from which obtain qualities
   * \param seq The Seqeunce from which obtain the characters
   *
   * \sa FullQuality(const QualifiedSequence& qualSeq)
   * \sa FullQuality(size_t n)
   */
  FullQuality(const Quality& qual, const Sequence& seq);
  /**
   * \brief Cosntruct a FullQuality from a QualifiedSequnce
   *
   * This constructor works as 
   * FullQuality(const Quality& qual,const Sequence& seq)
   * except that the Quality and the Sequence instances are taken
   * from the QualifiedSequence instance passed as parameter
   * \param qualSeq The QualifiedSequence from which retrieve Quality
   * and Sequence objects.
   *
   * \sa FullQuality(const Quality& qual, const Sequence& seq)
   * FullQuality(size_t n)
   */
  FullQuality(const QualifiedSequence& qualSeq);
  ~FullQuality();
  // ---------------------------------------------------------
  //                         OPERATORS                        
  // ---------------------------------------------------------

  /**
   * \brief Returns the probability vector for a given character
   *
   * This operator is used to retrieve the sequence of qualities
   * associated with the character `c`. The character is converted
   * to an index using the static method `getIndex` of the
   * template class C.
   * \param c The character for which we want the probabilities vector
   * \return a vector of proper size containing, for each poisiton of
   * the associated sequence, the probabilities of correctness of the
   * character c
   */
  double* operator[](char c);
private:
  // ---------------------------------------------------------
  //                     UTILITY FUNCTIONS                    
  // ---------------------------------------------------------
  double** initProbMatrix(size_t rows, size_t cols);
  void destroyProbMatrix(double** v, size_t rows, size_t cols);
};

/*************** CONSTRUCTORS AND DESTRUCTOR ***************/

template<class C>
FullQuality<C>::FullQuality(size_t n) {
  this->symbolCount = C::length();
  this->length = n;
  this->probVector = initProbMatrix(this->symbolCount, this->length);
}

template<class C>
FullQuality<C>::FullQuality(const Quality& qual, const Sequence& seq) {
  this->length = std::min(qual.length(), seq.getSequenceLength());
  this->symbolCount = C::length();
  this->probVector = initProbMatrix(this->symbolCount, this->length);
  // get the probabilities vector
  double* probs = qual.getProbabilities();
  // loop through all elements
  for (size_t pos = 0; pos < this->length; ++pos) {
    // get the index for the actual symbol
    size_t symbolIndex = C::getIndex(seq.getBaseAt(pos));
    double remind_p = (probs[pos] ) / ((double)(this->symbolCount - 1));
    // loop through all symbols
    for (size_t sym = 0; sym < this->symbolCount; ++sym) {
      this->probVector[sym][pos] = remind_p;
    }
    this->probVector[symbolIndex][pos] = 1.0 - probs[pos];
  }
}

template<class C>
FullQuality<C>::FullQuality(const QualifiedSequence& qualSeq)
  : FullQuality(*qualSeq.getQuality(), *qualSeq.getSequence())
{
}

template<class C>
FullQuality<C>::~FullQuality() {
  destroyProbMatrix(this->probVector, this->symbolCount, this->length);
}

/************************* OPERATORS ************************/

// typical usage would be
//     double qual_A_0 = fullQuality['A'][0];

template<class C>
double* FullQuality<C>::operator[] (char c) {
  return this->probVector[C::getIndex(c)];
}

/********************* UTILITY FUNCTIONS ********************/

template<class C>
double** FullQuality<C>::initProbMatrix(size_t rows, size_t cols) {
  double** v = new double*[rows];
  for (size_t i = 0; i < rows; ++i) {
    v[i] = new double[cols];
  }
  return v;
}

template<class C>
void FullQuality<C>::destroyProbMatrix(double** v, size_t rows, size_t cols) {
  for (size_t i = 0; i < rows; ++i ) {
    delete[] v[i];
  }
  delete[] v;
}

/************************************************************/

#endif
