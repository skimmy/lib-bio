#ifndef REFERENCE_H
#define REFERENCE_H

#include <core/Sequence.h>
#include <core/KMer.hpp>

#include <memory>
#include <cstring>
#include <list>

/**
 * \brief This class represents a \e refernce sequence from the
 * genomic point of view.
 *
 * References are genomic sequence (and therefore this class is
 * of Sequence type) and are usually used to represent known DNA
 * and or RNA sequences (less frequently there are other type of
 * sequences like proteic. amino acis, ...).
 *
 * This class is used to define a DNA biological sequence (although
 * it may be used also in other cases).
 *
 * \sa Sequence
 */
class Reference : public Sequence {

public:

  // -------------------------- STATIC FACTORY METHODS -------------------------

  static unique_ptr<Reference> createFromString(const string& s);

  // ----------------------- CONSTRUCTORS AND DESTRUCTOR  ----------------------
  /**
   * \brief Creates an empty Reference object
   *
   * \sa Reference(const string& sequence)
   * \sa Reference(const char* sequence, size_t n)
   * \sa Reference(const Reference& other)
   */
  Reference();
  /**
   * \brief Creates a Reference object from astandard \c string.
   *
   * This constructor implicitally assumes that the reference is
   * stored as a string character, the conversion is perfomed
   * via the <tt>c_str()</tt> method of \c string .
   * The size of the reference is infered from the length of the
   * passed string.
   * 
   * \param sequence The string containing the reference
   * 
   * \sa Reference(const Reference& other)
   * \sa Reference(const char* sequence, size_t n)
   * \sa Reference()
   */
  Reference(const string& sequence);
  /**
   * \brief Createsa a Reference object from a standard C string
   *
   * This constructor uses the input \c char pointer to create
   * a new Reference object with the given lenght. Input pointer
   * should refere to an array with at least the indicated number
   * of characters otherwise unpredictable behavior may arise.
   *
   * \param sequence The pointer to the sequence
   * \param n The length of the referenceto be created
   *
   * \sa Reference(const string& sequence)
   * \sa Reference()
   * \sa Reference(const Reference& other)
   */
  Reference(const char* sequence, size_t n);
  /**
   * \brief Creates a copy of the given Reference object
   *
   * \param other The Reference object to be copied
   *
   * \sa Reference(const string& sequence)
   * \sa Reference(const char* sequence, size_t n)
   * \sa Reference()
   */
  Reference(const Reference& other);
  ~Reference();

  // ---------------------------- GET AND SET METHODS ----------------------------
  /**
   * \brief Returns the list of all k-mers of the stored seuqence.
   *
   * This method creates and returns a standard \c list object
   * containing all the k-mers (i.e. k long substrings) of the
   * stored reference. The list is templated with KMer class and,
   * therefore, programmer should use this (preferred) class to
   * manage k-mers
   *
   * The length of the returned substrings (i.e. \c k ) is controlled
   * by the passed parameter. If \f$ n \f$ is the length of the 
   * reference, than the returned list contains exactly \f$ n - k + \f$
   * elements. Note that for efficiency and semantic reasons, the
   * returned list is <b>not a set</b> and therefore duplicated
   * k-mers are include as well. On the other hand this could be
   * usefull to maintain the reference to the 'source' of the k-mer
   * since the method guarantees that i-th returned element is
   * exactly the i-th k-mer of the reference (i.e. is the k long
   * substrings starting at position i of the reference).
   *
   * \param k The length of the returned substring
   * \return A \c list standard object containing all k-mers as
   * KMer objects
   *
   * \sa KMers
   */
  list<KMer> getKMerList(size_t k);

  // ---------------------------------------------------------
  //                     CONVERSION METHODS
  // ---------------------------------------------------------
  /**
   * \brief Converts the interal sequence into a color space
   * sequence
   * 
   * This method should be called only if necessary because it
   * requires creation of two new \c string object (to be passed
   * to the ColorAlphabet::basesToColors() method) used as
   * temporary container. For long sequences this operation may
   * sensibly slow down program execution.
   *
   * \param primer The \e primer character to initialize deconding
   * \return A referenced copy of the acutal Reference object (after
   * its conversion to color space representation)
   *
   * \sa toBases()
   */
  Reference& toColors(char primer);
  /**
   * \brief Converts the interal sequence into a base space
   * sequence
   * 
   * This method should be called only if necessary because it
   * requires creation of two new \c string object (to be passed
   * to the ColorAlphabet::colorsToBases() method) used as
   * temporary container. For long sequences this operation may
   * sensibly slow down program execution.
   *
   * \param primer The \e primer character to initialize deconding
   * \return A referenced copy of the acutal Reference object (after
   * its conversion to base space representation)
   *
   * \sa toColors()
   */
  Reference& toBases(char primer);

  // ---------------------------------------------------------
  //                'SEQUENCE' OVERRIDE METHODS
  // ---------------------------------------------------------
  const void* getSequence() const;
  size_t getSequenceLength() const;
  size_t getElementSize() const;
  size_t getByteCount() const;
  char getBaseAt(size_t i) const;

  // IOSTREAM
  friend ostream& operator<< (ostream& os, const Reference& ref);

private:
  // -------------------- UTILITY METHODS --------------------
  void init(size_t n);  


  char* sequence;
  size_t length;
};

#endif
