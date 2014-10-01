#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <string.h>

#include <cstddef>

#include <boost/functional/hash.hpp>
using namespace boost;

/**
 * \brief This is an abstract class representing a generic sequence.
 *
 * A sequence is characterized by its length and the size of
 * a single element.
 */
class Sequence {
 public:
  /**
   * \brief Returnsthe raw representation of the sequence.
   * 
   * The returned sequence is a raw copy of the original one,
   * it is intended to be used in a very general way, so safe
   * conversione can be performed only to those type that are
   * supported and the list of such types depend on the specific
   * implementation.
   * The length of the returned sequence (in bytes) can be obtained
   * using getByteCount() while the size of a single
   * element can be obtained calling getByteCount()
   *
   * \return The raw representation of the sequence.
   *
   * \sa getByteCount()
   * \sa getElementSize()
   */
  virtual const void* getSequence() const = 0;
  /**
   * \brief Returns the number of elements contained in the
   * sequence
   *
   * \return The number of elements contained in the sequence
   */
  virtual size_t getSequenceLength() const = 0;
  /**
   * \brief Returns the size of a single element.
   *
   * The size (until otherwise stated) should be intended
   * as the number of bytes.
   *
   * \return The size of a single element of the set
   */
  virtual size_t getElementSize() const = 0;
  /**
   * \brief Returns the size in byte of the sequence.
   *
   * \return The size in bytes of the sequence
   *
   * \sa getSequence()
   */
  virtual size_t getByteCount() const = 0;
  /**
   * \brief Returns the element at the given position
   *
   * \param i The position of the element to be returned
   * \return The element at the given position
   */
  virtual char getBaseAt(size_t i) const = 0;
};

/**
 * \brief This structure represents a unary function to convert a
 * Sequence into a hash
 *
 * This structure is used in combination with the
 * \e boost library to perform hasing on Sequence objects.
 */
struct sequence_hash : std::unary_function<Sequence, size_t> {
  /**
   * \brief Allows programmer (and the boost library) to perform
   * a call <tt>sequenceHash(sequnce)</tt> to obtain an hash value
   * for the given sequence.
   *
   * \param sequence The Sequence to be hashed
   * \return The hash value for the given equence
   */
  size_t operator()(const Sequence& sequence) const {
    size_t seed = 0;
    size_t n = sequence.getByteCount();
    char* rawSeq = (char*)sequence.getSequence();
    for (size_t i = 0; i < n; i++ ) {
      boost::hash_combine(seed, rawSeq[i]);
    }
    return seed;
  }
};


/**
 * \brief This structure represents a binary function to compare
 * two Sequence objects.
 *
 * The comparison is performed in two steps:
 * 1. The two sequences must have the same byte count as returned
 * by tge Sequence::getByteCount() method
 * 2. The two raw sequences returned by the Sequence::getSequence()
 * method must be equal compared with \c memcmp library function.
 */
struct sequence_equality : std::binary_function<Sequence, Sequence, bool> {
  /**
   * \brief Allows programmer (and boost library) to test equality betwwen
   *  two sequences with <tt>sequenceEquality(s1,s2)</tt>
   *
   * \param s1 The first Sequence to compare
   * \param s2 The second Sequence to compare
   * \return \c true if and only if the two input Sequence objects are equal
   *
   */
  bool operator()(const Sequence& s1, const Sequence& s2) const {
    return ( (s1.getByteCount() == s2.getByteCount()) && 
	     ( memcmp(s1.getSequence(), s2.getSequence(), s1.getByteCount()) == 0 ) );
  }
};



#endif
