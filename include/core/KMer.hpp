#ifndef K_MER_H
#define K_MER_H

#include <core/Sequence.h>

#include <cstdint>
#include <string>
#include <iostream>
using namespace std;

/**
 * \brief This class represents a \c k long sequence of characters
 *
 * This class is intended to be the representation of a \c k long
 * sequence (knwon as <em>k-mer</em>. The class provides basic facility
 * to assign and compare KMer onjects along with those methods inherited
 * from the Sequence interface
 *
 * \sa Sequence
 * \sa Read
 */
class KMer : public Sequence {
private:
  size_t k;
  string sequence;
public:
  // ---------------------------------------------------------
  //                CONSTRUCTOS AND DESTRUCTOR
  // ---------------------------------------------------------
  /**
   * \brief Creates an empty KMer with the given \c k
   *
   * \sa KMer(const string& seq)
   * \sa KMer(const char* seq, size_t k)
   * \sa KMer(const KMer& other)
   */
  KMer(size_t k);
  /**
   * \brief Creates a KMer from the given string.
   *
   * Both the sequence and \c k are determined by the input
   * string.
   *
   * \param seq The input sequence
   *
   * \sa KMer(const KMer& other);
   * \sa KMer(const char* seq, size_t k)
   * \sa KMer(size_t k)
   */
  KMer(const string& seq);
  /**
   * \brief Creates a KMer object from a standard C string.
   *
   * The passed reference should point to an array with (at
   * least) k elements otherwise unpredictable behavior may
   * arise. The reference should point to a greater array
   * but this constructor will simply ignore all characters
   * coming from the k-th position on.
   *
   * \param seq The sequence containing the k-mer
   * \param k The length of the k-mer
   *
   * \sa KMer(const string& seq)
   * \sa KMer(size_t k)
   * \sa KMer(const KMer& other);
   */
  KMer(const char* seq, size_t k);
  /**
   * \brief Creates a copy of the input KMer.
   *
   * \sa KMer(const string& seq)
   * \sa KMer(const char* seq, size_t k)
   * \sa KMer(size_t k)
   */
  KMer(const KMer& other);

  // ---------------------------------------------------------
  //                    GET AND SET METHODS
  // ---------------------------------------------------------
  /**
   * \brief 
   */
  string toString() const;

  // ---------------------------------------------------------
  //                         OPERATORS
  // ---------------------------------------------------------
  /**
   * \brief Writes the KMer on an output stream
   *
   * \param os The output stream
   * \param kmer The KMer to written
   * \return The output stream
   */
  friend ostream& operator<< ( ostream& os, const KMer& kmer);
  /**
   * \brief Assign operator
   */
  KMer* operator=(const KMer& other);
  /**
   * \brief Equality operator.
   *
   * The comparison between two KMer objects is performed by
   * comparing both the size \c k and the actual sequence
   *
   */
  bool operator==(const KMer& other);

  // ---------------------------------------------------------
  //                'SEQUENCE' OVERRIDE METHODS
  // ---------------------------------------------------------
  const void* getSequence() const;
  size_t getSequenceLength() const;
  size_t getElementSize() const;
  size_t getByteCount() const;
  char getBaseAt(size_t i) const;
};

#endif
