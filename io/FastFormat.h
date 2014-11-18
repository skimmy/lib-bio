#ifndef FASTFORMAT_H
#define FASTFORMAT_H

#include "Format.h"
#include "../sequence/Reference.hpp"

/**
 * \brief This class represents a <em>fast format</em> file.
 *
 * The class is should be used to load fast files containing a single
 * sequence (e.g. genomic and chromosomic sequences). Although it can
 * be used with any content.
 *
 * The parsing rules are very simple.
 * - A line of the file starting with the character \c > is
 * interpreted as an header line and the content is appended to
 * the header string (along with the newline).
 * - All other lines are interpreted as sequence lines and the
 * content of the line is appended to the sequence string (the
 * new line character is discarded).
 *
 * \sa Format
 */

class FastFormat : public Format {
 private:
  string sequence;
  string header;

 public:
  // ---------------------------------------------------------
  //                       CONSTRUCTORS
  // ---------------------------------------------------------
  /**
   * \brief Creates an empty FastFormat object
   *
   * The only operations perfromed by the default constructors are
   * those required by Format::Format() constructor requiring a name,
   * <em>FAST</em> is the ssigned
   *
   * \sa Fast::Format(const string &name)
   */
  FastFormat();

  
  /**
   * \brief Creates a FastFormat object by loading the fast file
   * with path passed as parameter.
   *
   * This is a shortcut for the invocation of default constructor
   * followed by a call to loadFromFile method.
   * 
   */
  FastFormat(const std::string& filePath);

  // ---------------------------------------------------------
  //                     CONVERSION METHODS
  // ---------------------------------------------------------
  /**
   * \brief Converts the FastFormat into a Reference objct
   *
   * \return A Reference object with the same content as the
   * actual FastFormat instance
   *
   * \sa FastFormat::operator Reference() const
   * \sa Reference
   */
  Reference toReference() const;
  
  // ---------------------------------------------------------
  //                         OPERATORS
  // ---------------------------------------------------------
  /**
   * \brief Converts the FastFormat into a Reference objct
   *
   * \return A Reference object with the same content as the
   * actual FastFormat instance
   *
   * \sa toReference()
   * \sa Reference
   */
  operator Reference() const;
  // ---------------------------------------------------------
  //                 'FORMAT' METHODS OVERRIDE
  // ---------------------------------------------------------
  string loadFromFile(const string &fileName);  
  string getSequence() const;
  string getHeader() const;
};

#endif
