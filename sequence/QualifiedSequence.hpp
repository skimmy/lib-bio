#ifndef QUALIFIED_READ_H
#define QUALIFIED_READ_H

#include "Sequence.h"
#include "../quality/Quality.hpp"

/**
 * \brief This class represents a sequence with quality 
 * values associated.
 *
 * The sequence is internally represented with a Sequence
 * pointer and the quality by a Quality pointer.
 * For efficiency reasons internal object are not allocated
 * and they are passed either by the constructor or using
 * the set methods. For this reasons all references should
 * remain valid during the entire scope of the object.
 *
 * \sa Quality
 * \sa Sequence
 */
class QualifiedSequence {
private:
  Quality* quality;
  Sequence* sequence;
public:
  // ---------------------------------------------------------
  //                CONSTRUCTORS AND DESTRUCTOR
  // ---------------------------------------------------------
  /**
   * \brief Initializes the QualifiedSequence with the two
   * pointer passes as paramaters
   *
   * \param seq The Sequence pointer
   * \param qual The Quality pointer
   */
  QualifiedSequence(Sequence* seq, Quality* qual);
  ~QualifiedSequence();
  // ---------------------------------------------------------
  //                    GET AND SET METHODS
  // ---------------------------------------------------------
  /**
   * \brief Returns the stored Quality pointer
   *
   * \return The stored Quality pointer
   *
   * \sa setQuality()
   * \sa getSequence()
   */
  Quality* getQuality() const;
  /**
   * \brief Returns the stored Sequence pointer
   *
   * \return The stored Sequence pointer
   * 
   * \sa setSequence()
   * \sa getQuality()
   */
  Sequence* getSequence() const;
  /**
   * \brief Sets the stored Quality pointer
   *
   * \param qual The \e new Quality pointer
   *
   * \sa setSequence()
   * \sa getQuality()
   */
  void setQuality(Quality* qual);
  /**
   * \brief Sets the stored Sequence pointer
   *
   * \param seq The \e new Sequence pointer
   *
   * \sa setQuality()
   * \sa getSequence()
   */
  void setSequence(Sequence* seq);
};

#endif
