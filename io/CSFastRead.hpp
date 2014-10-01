#ifndef CSFAST_READ_H
#define CSFAST_READ_H

#include "../sequence/Read.h"
#include "../sequence/KMer.hpp"
#include "../sequence/ColorAlphabet.hpp"
#include "../sequence/FullyQualifiedSequence.hpp"

#include "../quality/PhredQuality.hpp"


#include <iostream>
#include <list>
using namespace std;

/**
 * \brief This class represents a Read stored in the colors space format
 *
 * In <em>color space</em> sequence and quality are stored in two
 * different files, for this reason the preferred method to load
 * a read from input streams is the loadBasesAndQualitiesFromFiles()
 * method since it can handle two streams (rather the a single one
 * like the operator>>() does). To maintain symmetry on the stream
 * operators the operator<<() <b>doesn't write quality information</b>
 * whihc should be written in a different file to maintain the color
 * space standard.
 *
 * Recall, moreover, that reads in color space need a \e primer base
 * to be succesfully translated into base sequencec and, therefore,
 * methods to get and set such primer are also provided by the class.
 *
 * Although this class inherites get and set base methods from its
 * superclass Read, the sequnce that are actually stored in the
 * <em>base string</em> contains colors character (and not bases)
 * as also pointed out in the Read class description.
 *
 * \sa Read
 * 
 */
class CSFastRead : public Read {
private:
  string qualityFileName;
  char primer;
public:
  // ---------------------------------------------------------
  //                      CONSTRUCTORS
  // ---------------------------------------------------------

  /**
   * \brief Creates an empty CSFastRead
   *
   * This constructor creates and emtpy color space read invoking
   * the Read::Read() superclass constructor.
   *
   * \sa CSFastRead(const CSFastRead& other)
   * \sa Read::Read()
   *
   */
  CSFastRead();
  /**
   * \brief Creates a CSFastRead object from an existing one
   *
   * The copy constructor creates an exact copy of the given
   * CSFastRead. Since all internal variables are copied into
   * new one, this constructor should be avoided when massive
   * reads set must be managed both for space and time efficiency.
   *
   * \param other The CSFastRead to be copied
   *
   * \sa CSFastRead()
   */
  CSFastRead(const CSFastRead& other);


  // ---------------------------------------------------------
  //                    GET AND SET METHODS
  // ---------------------------------------------------------
  /**
   * \brief Returns the primer base actually stored in the color
   * space read
   *
   * \return The prime base
   *
   * \sa setPrimer()
   */
  char getPrimer() const;
  /**
   * \brief Changes the prime base of the CSFastRead
   *
   * \param p The new primer base
   *
   * \sa getPrimer()
   */
  void setPrimer(char p);

  /**
   * \brief Returns a standard \c list of FullyQualifiedSequence
   * one for each k-mer of the read.
   *
   * The list is created starting from the stored seuqence and the
   * stored quality values. The class requires that the sequence is
   * represented into color space (since ColorAlphabet class will
   * be used) and that qualities a represented as phred values.
   * This assumptions are reasonable due to the standard format
   * of \c csfast files.
   * The returned list contains \f$ n-k+1 \f$ elements (where \c n
   * is the size of the read and \c k ise the size of k-mers), each
   * element is a FullyQualifiedSequence templated on ColorAlphabet.
   * Recall that the references passed to constructors of the
   * FullyQualifiedSequence class require that <b>pointer to the
   * Sequence and FullQuality objects remain valid through all the
   * class scope</b> and therefore they are dynamically allocated
   * and <b>the programmer should take care of freeing them</b>
   * when they are no more needed.
   *
   * \param k The size of the k-mers
   * \return A list of FullyQualifiedSequence for the k-mers
   *
   * \sa Read::getKMerList()
   * \sa FullQuality
   * \sa Sequence
   * \sa FullyQualifiedSequence
   * \sa ColorAlphabet
   * \sa PhredQuality
   *
   */
  list< FullyQualifiedSequence< ColorAlphabet > > getFullyQualifiedKMerList(size_t k) const;

  // ---------------------------------------------------------
  //                  '>>' AND '<<' OPERATORS
  // ---------------------------------------------------------
  /**
   * \brief Loads the header and color sequence from an input
   * stream.
   *
   * Since color sequence and quality values are stored in
   * different files this operator (which can take only one
   * input stream as parameter) loads only the color sequnce
   * and the header of the read. To load both color sequnce
   * and quality values use loadBasesAndQualitiesFromFiles().
   *
   * \param is The input stream
   * \param read The CSFastRead where the loaded read should
   * be stored
   * \return The input stream used to load the read
   *
   * \sa loadBasesAndQualitiesFromFiles()
   * \sa operator<<()
   */
  friend istream& operator>>(istream& is, CSFastRead& read);
  /**
   * \brief Writes on an ouptut stream the color sequence
   * and the header of the read
   *
   * To maintain \e consistency with the operator>>(), this
   * operator writes only the header and the color sequence
   * (along with the stored primer).
   *
   * \param os The output stream
   * \param read The read to write on the stream
   * \return The ouptu stream
   *
   * \sa operator>>()
   */
  friend ostream& operator<<(ostream& os, CSFastRead& read);

  /**
   * \brief Loads color sequence, header and quality values
   * from two input streams.
   *
   * This should be the preffered method to load a read from
   * input files. This methods assumes that both streams are
   * \e syncrhonized in the sense that they point to the same
   * read (i.e. the read with the same header). If this constrain
   * is violated, the method doesn't load the read but, doesn't
   * terminate the execution. Repeatdly calling this method
   * fusing the same streams (poitning to well formed color
   * space files) guarantees the success of load operations.
   * 
   * \param bases_is The input stream pointing to the sequence file
   * \param qual_is The input stream pointing to the qualiy file
   *
   * \sa operator<<()
   * \sa operator>>()
   * 
   */
  void loadBasesAndQualitiesFromFiles(istream& bases_is, istream& qual_is);
  // TODO: provide a writeBasesAndQualitiesToFiles() method
};

#endif
