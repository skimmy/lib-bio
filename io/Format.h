#ifndef FORMAT_H
#define FORMAT_H

#include <string>
using namespace std;

/**
 * \brief This class represents a file format for storing sequences.
 *
 * The class is designed to be used within the bioinformatics field
 * and, therefore, it contains all primitives to access all
 * informations usulally stored in such files. These informations
 * are usually of two kinds. The first is the actual sequence and
 * the secondo is a (possibly empty) header describing the sequence
 * itself.
 */
class Format {
 protected:
  /**
   * \brief The name of the format
   */
  string name;
 public:
  /**
   * \brief Initialize the Format with a name
   *
   * This is the only constructor provided the class and it simply
   * initializes the name of the format.
   */
  Format(const string &name);
  /**
   * \brief Loads from file
   *
   * This is virtual function beacause the actual content of the
   * file depends on the specific format
   */
  virtual string loadFromFile(const string &fileName) = 0;  
  /**
   * \brief Returns the sequence
   *
   * \return The sequence
   */
  virtual string getSequence() const = 0;

  /**
   * \brief Returns the header
   *
   * \return The header
   */
  virtual string getHeader() const = 0;
};

#endif
