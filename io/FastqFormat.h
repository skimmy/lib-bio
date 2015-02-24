#ifndef FASTQ_FORMAT_H
#define FASTQ_FORMAT_H

#include "../io.h"

#include <string>
#include <iostream>
#include <fstream>
using namespace std;

/**
 * \brief This class represents a <em>Fastq format</em> file.
 *
 * Fastq format maintains in a single file both sequence and
 * quality information. This class has been designed for work
 * best with reads fastq files (i.e. files containing multiple
 * sequence each with an header, a base sequence and a quality
 * tring).
 *
 * The class is designed to be used in combination with FastqRead
 * file. Using the getNextRead() method programmers may access
 * the next read stored in the file. Nevertheless it is possible
 * to load the entire file into the internal string using the
 * loadFromFile() method, this approach is however discouraged
 * for huge because program may rapidly run out of memory (and
 * eventually crashing or slowing impressively).
 *
 * \sa Format
 * \sa FastqRead
 */
class FastqFormat : public Format {
  // private fields
 private:
  string sequence;
  string header;
  ifstream inFile;
  FastqRead* nextRead;


 public: 
  // ---------------------------------------------------------
  //                CONSTRUCTORS AND DESTRUCTOR
  // ---------------------------------------------------------
  /**
   * \brief Creates an empty FastqFormat object.
   *
   * The default constructor simply initialize its name with
   * the string \e FASTQ (as prescribed by the Format class),
   * while all other internal properties are set to empty.
   *
   * \sa Format::Format()
   */
  FastqFormat();
  ~FastqFormat();
 
  // ---------------------------------------------------------
  //                       LOAD METHODS
  // ---------------------------------------------------------
  /**
   * \brief Opens a file containing fastq sequences
   *
   * The method doesn't perform any check on the availability
   * and on the content of the file, it simply returns a \c bool
   * indicating whether or not opening an input stream on the
   * indicated file has been succesfull
   *
   * \param fileName The full path of the fastq file
   * \return \c true if opening an input stream on the file
   * succeded
   *
   * \sa loadFromFile()
   */
  bool openFile(const string& fileName);
  /**
   * \brief Check whether the file contains more reads
   *
   * Once an input file is openede with the openFile() method
   * the class is able to inspect whether or not there is more
   * to be read from the file. Note that finding more characters
   * on the input file doesn't mean that actually the characters
   * represent a read, this check is (for efficiency reasons)
   * not performed and malformed input files could arise
   * unpredictable behaviors.
   *
   * \return \c code if the file contains more characters
   *
   * \sa getNextRead()
   * \sa openFile(const string& fileName)
   */
  bool hasNextRead();
  /**
   * \brief Loads the next read from the input file (if any)
   *
   * This method must be used in combination with hasNextRead()
   * and should be called only after a succesfull call to openFile()
   * has been performed.
   * If the file still contains a read it is actually returned by
   * this method. If the file dosen't contain any character (i.e.
   * hasNextRead() would return false) anm empty read is returned.
   * If the input file contains more characters but these characters
   * don't represent a well formed fastq read, the behavior could
   * be unpredictable.
   *
   * \return The new read loaded from file (or an empty one if the
   * input file doesn't contain any character)
   *
   * \sa hasNextRead()
   * \sa openFile(const string& fileName)
   */
  FastqRead getNextRead();

  // ---------------------------------------------------------
  //                 'FORMAT' OVERRIDE METHODS
  // ---------------------------------------------------------
  string loadFromFile(const string &fileName);
  string getSequence() const;
  string getHeader() const;

 private:
  // ---------------------------------------------------------
  //                      UTILITY METHODS
  // ---------------------------------------------------------
  void loadNextRead();
};

#endif
