#ifndef COMPRESSED_SEQUENCE_H
#define COMPRESSED_SEQUENCE_H

#include "Sequence.h"

#include <cstdint>
#include <cstddef>
#include <string>
using namespace std;

/**
 * \brief This class represents a sequence of \em compressed elements.
 *
 * Elements are compressed in the sense that, given the size in bits
 * of single element, the minumum number of \b bytes are used to
 * store the entire sequence. If the number of bits is less than
 * one byte more than one element is stored per byte.
 *
 * The class is designed to be used with element sizes (in bits)
 * that are power of two: 1,2,4 and 8 are typycal values although
 * 8 bits elements are best represented using a sequence of bytes.
 *
 * Internally the sequence is stored as an array of bytes, if
 * \f$ n \f$ is the number of element in the sequence and \f$ m \f$
 * is the size (in bits) of a single element then the real size 
 * \f$ N \f$ of the internal sequence is
 * \f[
 * N = \lceil \frac{n}{8} \cdot m \rceil
 * \f]
 * This size can be obtained by calling the length() method and
 * is usually refered to as the <em>real size</em> of the sequence.
 */
class CompressedSequence : public Sequence {
 protected:
  /**
   * \brief The raw sequence
   */
  uint8_t* seq;
  /**
   * \brief The number of elements of the sequence
   */
  size_t n;
  /**
   * \brief The real size (in bytes) of the sequence
   */
  size_t realSize;
  /**
   * \brief The size (in bits) of the single element
   */
  size_t elSize;
  /**
   * \brief The mask used to access single element within 
   * a whole byte
   */
  uint8_t mask;
 public:
  // ---------------------------------------------------------
  //                CONSTRUCTORS AND DESTRUCTOR
  // ---------------------------------------------------------


  /**
   * \brief Default constructor, an zero length sequence of 
   * elements of zero bytes size is created.
   *
   * This constructor can be used in order to create an
   * instance which will be assigned <em>in the future</em>
   *
   * \sa CompressedSequence(size_t n, size_t elSize)
   * \sa CompressedSequence(const CompressedSequence& s)
   * 
   */
  CompressedSequence();

  /**
   * \brief Creates a new CompressedSequence with the given
   * number of elements and size of a single element.
   *
   * The number of elements in the sequence is the first
   * argument and it is mandatory. The size (in bits) of a single
   * element is the second parameter and can be omitted in
   * which case the default value is 2 bits
   * 
   * \param n The number of elements in the sequence
   * \param elSize The size (in bits) of a single element
   *
   * \sa CompressedSequence()
   * \sa CompressedSequence(const CompressedSequence& s)
   */
  CompressedSequence(size_t n, size_t elSize = 2);

  /**
   * \brief Copy constructor to initialize from another 
   * CompressedSequence instance.
   *
   * \param s The CompressedSequence to be copied
   *
   * \sa CompressedSequence()
   * \sa CompressedSequence(size_t n, size_t elSize)
   */
  CompressedSequence(const CompressedSequence& s);
  ~CompressedSequence();

  // ---------------------------------------------------------
  //                    GET AND SET METHODS
  // ---------------------------------------------------------

  /**
   * \brief Returns the <em>real size</em> of the sequence.
   *
   * The real size is defined as the number of bytes that the
   * sequence actually occupy in its internal representation.
   * The same result can be obtained calling the getByteCount(),
   * this redundancy is maintained to give more expressivity to
   * the class while obeying the hierarchy prescriptions (inherited
   * from the Sequence interface)
   * 
   * \return The real size of the sequence
   *
   * \sa getElementCount() const
   * \sa getSequenceLength() const
   * \sa getElementSize() const
   * \sa getByteCount() const
   */
  size_t length() const;

  /**
   * \brief Returns the number of element stored in the sequence.
   *
   * \return The number of elements stored
   *
   * \sa length() const
   * \sa getElementSize() const
   * \sa getSequenceLength() const
   * \sa getByteCount() const
   *
   */
  size_t getElementCount() const;

  /**
   * \brief Returns the i-th element of the sequence.
   *
   * The returned is a byte which contains in the less
   * significant bits the actual element. In other words
   * for an element size of 2 only the two less significant
   * bit are set to the proper value while the remaining 6
   * bits are set to zero.
   *
   * \param i The position from which we want the element
   * \return The element in the i-th position
   * 
   * \sa getBaseAt(size_t i)
   * \sa setElementAt(size_t i, uint8_t e)
   */
  uint8_t getElementAt(size_t i) const;

  /**
   * \brief Sets the i-th element of the sequence.
   *
   * The setting of the element take care of considering the
   * proper number of bits from the parameter.
   * \param i
   * \param e
   *
   * \sa getElementAt(size_t i)
   * \sa getBaseAt(size_t i)
 */
  void setElementAt(size_t i, uint8_t e);

  /**
   * \brief Returns the internal sequence in its \e raw format.
   *
   * The method returns the sequence as a pointer, the returned
   * reference points to <b>the same sequence stored internally</b>
   * and therefore any change using the returned pointer will
   * modify the internal sequence of the class.
   *
   * The only difference with the getSequence() method is that
   * getawSequence() returns a uint8_t pointer (that is also
   * the way the sequence is stored internally)
   * 
   * \return A pointer to the internal sequence
   *
   * \sa getSequence() const
   */
  const uint8_t* getRawSequence() const;

  // ---------------------------------------------------------
  //                      MODIFY METHODS
  // ---------------------------------------------------------

  /**
   * \brief Append another CompressedSequence to the actual one.
   *
   * The methods appends to the internal sequence the sequence
   * stored in the parameter. It assumes that both sequence are
   * \e homogeneous in the sense that have the same size of the
   * single element (namely the size stored on the actual
   * CompressedSequence)
   *
   * The method also returns a referenced copy of the actual
   * CompressedSequence after the append operation has been
   * performed. This is done in order to allow chain calls to
   * the same Compressed Sequence like
   * \code
   * s1.append(s2).append(s3);
   * \endcode
   *
   * \param other The CompressedSequence to be appended
   * \return The actual sequence after the append operation has
   * been performed
   */
  CompressedSequence& append(const CompressedSequence &other);

  // ---------------------------------------------------------
  //                        I/O METHODS
  // ---------------------------------------------------------

  /**
   * \brief Writes the sequence to a binary file
   *
   * In order to write the full information the output file
   * contains (at the beginning of the file) the binary formato
   * of the <em>real size</em> and of the <em>element count</em>.
   * After such information the raw sequence is written in
   * binary format.
   *
   * \param fileName The full path of the output file
   * 
   * \sa loadFromFile(const string& fileName)
   */
  void writeToFile(const string& fileName) const;

  /**
   * \brief Loads the sequence from a binary file
   *
   * The method assumes that the input binary file has the
   * format produced by the writeToFile() method. That is, 
   * at the beginning of the file it looks for the <em>real
   * size</em> and for the <em>element size</em> (as size_t
   * types) and afterward the raw sequence is loaded
   *
   * \param fileName The full path of the file to be loaded
   * \return a referenced copy of the sequence after the load
   * operation has been performed
   *
   * \sa writeToFile(const string& fileName) const
   */
  CompressedSequence& loadFromFile(const string& fileName);

  // 'Sequence' class override
  const void* getSequence() const;
  size_t getSequenceLength() const;
  size_t getElementSize() const;
  size_t getByteCount() const;
  char getBaseAt(size_t i) const;

 protected:
  // ---------------------------------------------------------
  //                      UTILITY METHODS
  // ---------------------------------------------------------

  /** 
   * \brief Initializes the sequence once all internal members
   * are set
   */
  void init();

  /**
   * \brief Resizes the internal sequence
   *
   * \param newSize The new size of the internal sequence
   */
  void resize(size_t newSize);

};

#endif
