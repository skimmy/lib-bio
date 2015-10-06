#ifndef CS_FAST_FORMAT_H
#define CS_FAST_FORMAT_H

#include "CSFastRead.hpp"
#include "Format.h"

#include <fstream>
using namespace std;

/**
 * \brief This class is an implementation of the Format class
 * for <em>Color Space</em> read input.
 *
 * Recall that color space reads are given in two different 
 * files: one containing the actual sequence of reads, and
 * another containing the quality values associated to the
 * calls of the read.
 *
 * The class requires the two input files (one for colors
 * and one for qualities) and allows the programmer to load
 * reads one at time using the getNextRead() method which
 * returns a CSFastRead object.
 *
 * \sa CSFastRead
 *
 */
class CSFastFormat : public Format {
 private:
  string sequence;
  string header;
  ifstream baseFile;
  ifstream qualFile;
  CSFastRead* nextRead;
public:
  // ---------------------------------------------------------
  //                 CONSTRUCTORS AND DESTRUCTOR
  // ---------------------------------------------------------

  /**
   * \brief This constructor associate the format with the two input files.
   * 
   * \param basesFileName The full path of the file containing the colors
   * \param qualFileName The full path of the file containint the qualities
   */
  CSFastFormat(const string& basesFileName, const string& qualFileName);
  ~CSFastFormat();

  // ---------------------------------------------------------
  //                      LOAD METHODS
  // ---------------------------------------------------------

  /**
   * \brief Check if there are still reads on the files.
   *
   * This method check whether or not there are more reads
   * to be loaded from the input files. For efficiency reson
   * the test is performed by simply checking whether or not
   * there are characters on the input streams and, therefore,
   * malformed files could result in a unpredictable result.
   *
   * \return true if there are still reads to be loaded from
   * the files, false otherwise
   *
   * \sa getNextRead()
   */
  bool hasNextRead();

  /**
   * \brief Loads a new color space read from the input files.
   *
   * Since no check is performed on the validity of the input
   * files, in case of malformed input, the result is
   * unpredictable. If there aren't more reads on file the
   * method returns an empty CSFastRead object.
   *
   * \return A new CSFastRead object from the input files
   * (or an empty one if no more reads are contained in the files).
   *
   * \sa hasNextRead()
   * \sa CSFastRead
   *
   */
  CSFastRead getNextRead();

  // Format methods override
  string loadFromFile(const string &fileName);  
  string getSequence() const;
  string getHeader() const;

private:
  // utility methods
  void loadNextRead();
};

#endif
