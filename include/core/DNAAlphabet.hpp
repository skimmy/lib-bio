#ifndef DNA_ALPHABET_H
#define DNA_ALPHABET_H

#include <string>
#include <map>
using namespace std;

/**
 * \brief This class contains static methods to represent a
 * <em>DNA Alphabet</em>
 *
 * This class is used to represent tha DNAAlphabet made of the
 * fourc character A, C, G and T. The class intentionally
 * doesn't implement a wider set (like DNA alphabet allowing
 * gaps or a more general IUPAC encoding alphabet) it is
 * however sufficient to represent \e ungapped sequences of
 * DNA bases.
 *
 * An alphabet is a set of characters (or any comparable object
 * in a more wide sense). Some classes (like FullQuality and
 * FullyQualifiedSequence) must be templated using alphabet
 * classes. An alphabet class should supply the basic methods
 * to convert characters into indexes and vice-versa.
 *
 * \sa ColorAlphabet
 * \sa FullQuality
 * \sa FullyQualifiedSequence
 */
class DNAAlphabet {
private:
  static const size_t LENGTH = 4;
public:
  /**
   * \brief Returns the length of the alphabet.
   *
   * \return The length of the alphabet
   */
  static size_t length();
  /**
   * \brief Converts an index into a character of the
   * alphabet
   *
   * This methods converts an index between 0 and
   * n-1 (where n is the size of the alphabet).
   * It can return both upper and lower case characters since in
   * DNA alphabet character case should be considered irrilevant
   *
   * \param  i The index to converto into character
   * \param upper When true the returned character is upper
   * cased, when false the returned character is lower cased
   * (default value is true)
   * \return The character corresponding to the input index
   *
   * \sa getIndex()
   * \sa length()
   *
   */
  static char getChar(size_t i, bool upper = true);
  /**
   * \brief Converts a character of the alphabet into an index.
   *
   * The character given as input could be upper or lower case
   * and should be between 0 and n-1 (where
   * n is the size of the alphabet) otherwise an error
   * would occur.
   * 
   * \param c The character to be converted into index
   * \return The index associated with the input character
   * 
   * \sa getChar()
   * \sa length()
   */
  static size_t getIndex(char c);
private:
  static map<char,size_t> initCharMap();
};

#endif
