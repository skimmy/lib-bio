#ifndef READ_H
#define READ_H

#include "DNACompressedSymbol.h"
#include "Sequence.h"
#include "KMer.hpp"

#include <string>
#include <list>


/**
 * \brief This class represents a read in the bioinformatics context.
 *
 * A read is made of three elements
 *  -# a sequence of bases,
 *  -# an optional header and
 *  -# a sequence of optional quality values.
 * There are no restrictions on the content of <em>bases, quality
 * and header</em> all of which could eventually be empty. The
 * length of the Read is determined by the number of characters
 * composing the `bases` field of the Read, the length so obtained
 * is also used to infer the length of the quality string (i.e.
 * the number of quality values in the string rather than the
 * length in character of the string). For qualitry sequences
 * containing less than such length the class behavior could
 * be unpredictable, while for sequences longer than the actual
 * length remaining values usually are ignored (although they
 * are still stored and returned whenever the getQualities()
 * method is called.
 *
 * Note that the term \em bases for the actual content of the
 * sequence is used for convinence, but read can be made of
 * colors or amino acid as well without inficing the semantic
 * of the class.
 */
class Read : public Sequence {
 protected:
  /**
   * \brief The header of the read
   */
  std::string header;
  /**
   * \brief The bases of the read
   */
  std::string bases;
  /**
   * The string containing quality values
   */
  std::string qualities;
 public:
  /**
   * \brief Creates an empty Read
   *
   * An empty read has emtpy strings for bases, qualities and
   * header.
   */

  Read();

  // ---------------------------------------------------------
  //                    SET AND GET METHODS
  // ---------------------------------------------------------
  /**
   * \brief Sets the header of the Read
   *
   * \param header The \em new header
   */
  void setHeader(const std::string& header);
  /**
   * \brief Sets the base sequence of the Read
   *
   * \param bases The new sequence
   */
  void setBases(const std::string& bases);
  /**
   * \brief Sets the qualities string of the Read
   *
   * \param qualities The new quality string
   */
  void setQualities(const std::string& qualities);

  /**
   * \brief Returns the header of the Read
   *
   * \return The stored header
   */
  std::string getHeader() const;
  /**
   * \brief Returns the sequence of bases of the Read
   *
   * \return The stored base sequence
   */
  std::string getBases() const;
  /**
   * \brief Returns the quality string of the Read
   *
   * \return The stored quality string
   */
  std::string getQualities() const;
  /**
   * \brief Returns the length (in bases) of the string
   *
   * \return The length of the stored base sequence
   */
  size_t length() const;

  /**
   * \brief Returns a standard list made of all k-mers.
   *   
   * K-mers all the substrings  of length k that are obtained
   * the sequence of bases of the Read. The method returns a
   * standard list object containing all the \f$ n - k + 1 \f$
   * k-mers of the read (where \f$ n \f$ is the number of bases
   * contained in the Read)
   *
   * \param k the size of the k-mer
   * \return A standard list containing all k-mers.
   *
   */
  std::list<KMer> getKMerList(size_t k) const;

  // ---------------------------------------------------------
  //                     MODIFY OPERATION
  // ---------------------------------------------------------
  /**
   * \brief Trims the read to the indicated size.
   *
   * This operation simply trims the base sequence
   * while leave unaltered header and qualities. While header is
   * supposed not to change when trimmed, qualities should change
   * accordingly, this operation is left as an operation performed
   * by the subclasses since it's dependent on the format specified
   * for the Read.
   *
   * \param n the new size of the read
   * \return this object after trimming has took place
   */
  Read& trim(size_t n);

  // ---------------------------------------------------------
  //                 'SEQUENCE' CLASS OVERRIDE
  // ---------------------------------------------------------
  const void* getSequence() const;
  size_t getSequenceLength() const;
  size_t getElementSize() const;
  size_t getByteCount() const;
  char getBaseAt(size_t i) const;
};

#endif
