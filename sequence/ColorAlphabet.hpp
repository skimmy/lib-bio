#ifndef COLOR_ALPHABET_H
#define COLOR_ALPHABET_H

#include <string>
#include <map>
using namespace std;

/**
 * \brief This class represents the <em>color space</em> alphabet.
 * 
 * The color space alphabet made of the four characters 0,1,2
 * and 3. The class doesn't consdier the . (dot) character used
 * in the color space reds to represent a gap.
 */ 
class ColorAlphabet {
private:
  static const size_t LENGTH = 4;
public:
  /**
   * \brief Returns the length of the color alphabet
   *
   * \return The length of the color alphabet
   */
  static size_t length();
  /**
   * \brief Converts an index into a character
   *
   * \param i The index to be converted
   * \return The char corresponding to the input index
   */
  static char getChar(size_t i);
  /**
   * \brief Converts a character to an index
   *
   * \param c the char to be covnerted
   * \return The index corresponding to the input character
   */
  static size_t getIndex(char c);
  /**
   * \brief Converts a color string into a base string
   *
   * This method takes a colors string and a primer character
   * to transform the sequence into a <em>base space</em> string.
   * 
   * \param colors The color string to be converted
   * \param bases The base output string 
   * \param primer The \e primer character for sequence decoding
   */
  static void colorsToBases(const string& colors, string& bases, char primer);
  /**
   * \brief Converts a base string into a color string
   *
   * This method takes a bases string and a primer character
   * to transform the sequence into a <em>color space</em> string.
   * 
   * \param bases The base string to be converted
   * \param colors The color output string 
   * \param primer The \e primer character for sequence decoding
   */
  static void basesToColors(const string& bases, string& colors, char primer);
private:
  static map<char,size_t> initCharMap();
};

#endif
