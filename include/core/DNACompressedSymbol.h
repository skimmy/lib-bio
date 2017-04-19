#ifndef DNA_COMPRESSED_SYMBOL
#define DNA_COMPRESSED_SYMBOL

#include <cstdint>
#include <cstddef>

/**
 * \brief This class represents a DNA symbol (i.e. character)
 *  in a compressed form.
 *
 * DNA symbols as defined by the IUPAC encoding are characters
 * representing any possible subset of the four DNA bases (A,
 * C, and T). This class allows programmer to manage such symbols
 * in a <em>compressed form</em> so that they occupy the less
 * possible space (in terms of number of bits). Namely for
 * a four character DNA alphabet we have 16 possible subsets and
 * therefore the minimum size of a single is 4 bits.
 *
 * This class is intended to be used in combination with the
 * CompressedSequence and CompressedReadSet classes sihce 
 * the ouput represantion is stored in the four less significant
 * bit of a uint8_t type (i.e. a byte).
 *
 * This class is also a perfect candidated (with slight
 * modifications to become an alphabet class for compressed
 * IUPAC alphabet as an extension of DNAAlphabet)
 *
 * \sa CompressedSequence
 * \sa CompressedReadSet
 * \sa DNAAlphabet
 * 
 */
class DNACompressedSymbol {
 private:
  uint8_t symbol;
 public:
  // ---------------------------------------------------------
  //                      CONSTRUCTORS
  // ---------------------------------------------------------
  /**
   * \brief This construcotr creates a new symbol from the indicated
   * character (which a \e gap character when omitted)
   *
   * \param c The character of the IUPAC encoding to be stored
   * in the DNACompressedSymbol (default is '.' meaning a gap)
   */
  DNACompressedSymbol(char c = '-');

  // ---------------------------------------------------------
  //                         OPERATORS
  // ---------------------------------------------------------
  /**
   * \brief Cast operator to convert DNACompressedSymbol into char
   * 
   * \return The conversion to a char
   *
   * \sa operator uint8_t()
   * \sa NumberToIupac()
   */
  operator char() const;
  /**
   * \brief Cast operator to convert DNACompressedSymbol into
   * a uint8_t type
   *
   * \return The conversion to uint8_t
   *
   * \sa operator char()
   * \sa IupacToNumber()
   */
  operator uint8_t() const;
  /**
   * \brief Assign operator
   */
  DNACompressedSymbol& operator = (const char c);

  // ---------------------------------------------------------
  //                       STATIC METHODS
  // ---------------------------------------------------------
  /**
   * \brief Converts a compressed representation into a character
   *
   * This is a convinience method to convert the compressed
   * representation into the corresponding character.
   * This method should be used when only the conversion is
   * needed and not an instance of the class
   *
   * \param n The compressed representation to be converted
   * \return The character representation
   *
   * \sa IupacToNumber()
   */
  static  char NumberToIupac(uint8_t n);
  /**
   * \brief Converts a character into its compressed form
   *
   * This is a convinience method to convert a character
   * into its compressed representation.
   * This method should be used when only the conversion is
   * needed and not an instance of the class
   *
   * \param c The character to be converted
   * \return The compressed representation of the input character
   *
   * \sa NumberToIupac()
   */
  static  uint8_t IupacToNumber(char c);
};

#endif
