#ifndef FASTQ_READ_H
#define FASTQ_READ_H

#include <core/Read.hpp>
#include <core/Quality.hpp>

#include <iostream>
#include <vector>

/**
 * \brief This class represents a read with a <em>fastq</em> format.
 * 
 * The class provides basic facilities to manage <em>fastq format</em>
 * reads. It is designed to be used in combination with the FastqFormat
 * class.
 *
 * Basic parsing methods are provided by means of the operator>>()
 * methos. The parsing performed by this class is very simple, recalling
 * the general fastq read format
 * \code
 * @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
 * GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC
 * +
 * SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
 * IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9IC 
 * \endcode
 * The rules are the following
 * - If the line starts with a \c \@ character it is interpreted as
 * the header (only single line header is allowed)
 * -  the line immedialty following the header is interpreted as
 * the base sequence
 * - The line between sequence and quality string (the \c + line
 * in our example) is discarded
 * - The last line is interpreted as a the quality string.
 * Although these parsing rules may fail when loading long sequences
 * (i.e. sequences spanning through multiple lines) the class works
 * fine with short sequences (i.e. those produced by <em>Next
 * Generation %Sequence</em> machines).
 *
 * A final remark regards the presence of quality strings which are expected to
 * be presented because prescribed by the `fastq` format. Such quality scores are
 * supposed to be represented on the phred scale. An attempt to determine whether
 * _sanger_ or _illumina_ encoding is used is automatically performed, if this
 * fails (and an expliciti indication is not given) sanger format is assumed.
 *
 * \sa Read
 * \sa FastqFormat
 * \sa Format
 * \sa CSFastRead
*/
class FastqRead : public Read {
 public:
  // ---------------------------------------------------------
  //                      CONSTRUCTORS
  // ---------------------------------------------------------
  /**
   * \brief Construct an empty FastqRead object
   *
   * \sa FastqRead(const FastqRead& other)
   */
  FastqRead();
  /**
   * \brief Create a copy of an existing FastqRead
   *
   * \sa FastqRead()
   */
  FastqRead(const FastqRead& other);

  ~FastqRead();

  // ---------------------------------------------------------
  //                 '>>' AND '<<' OPERATORS
  // ---------------------------------------------------------
  /**
   * \brief Loads a FastqRead from an input stream
   *
   * This operator allows the programmer to use an input
   * stream pointing to the beginning of a well formed
   * read entry (i.e. as described in the class description)
   * to load the content of the file
   * 
   * This operator is used by the FastqFormat class to
   * load a FastqRead from a file
   *
   * \param is The input stream
   * \param read The FastqRead object to be loaded
   * \return The input stream
   *
   * \sa operator<<()
   */
  friend std::istream& operator>>(std::istream& is, FastqRead& read);
  /**
   * \brief Writes the FastqRead on an output stream
   *
   * This operator writes the FastqRead object into an
   * output stream. The format of the output is coherent
   * with that described in the class description
   * 
   * \sa operator>>()
   */
  friend std::ostream& operator<<(std::ostream& os, const FastqRead& read);

  // ---------------------------------------------------------
  //                        STATIC METHODS
  // ---------------------------------------------------------
  static std::vector< size_t > splitReads(const std::string& filePath, size_t T);

  // ---------------------------------------------------------
  //                        QUALITY METHODS
  // ---------------------------------------------------------
  
  bool hasQualities() const { return (this->qualities.size() > 0); }

  /**
   * \brief Trys to decode quality string (if any) based on its content,
   */
  void autoDecodeQualities();

  bool hasDecodedQuality() const { return (this->quality != NULL); }

  Quality* getQuality() { return this->quality; }

 private:
  Quality* quality;
  
};

#endif
