#ifndef FULLY_QUALIFIED_SEQUENCE_H
#define FULLY_QUALIFIED_SEQUENCE_H

#include "FullQuality.hpp"

/**
 * \brief This class is actually a container for a Sequence
 * instance and FullQuality instance.
 *
 * This class contains an instance of Sequence and an instance
 * of FullQuality. The class should be templated with the same
 * type as the FullQuality instance. The class doesn't create
 * a copy of the passed reference, therefore the programmer must
 * ensure that passed references remain valid as long as they
 * will be needed (unpredictable behavior would otherwise arise).
 * Note that, moreover, no delete operation is performed by the
 * class' destructor.
 */
template<class C>
class FullyQualifiedSequence {
private:
  Sequence* sequence;
  FullQuality<C>* quality;
public:
  // ---------------------------------------------------------
  //                CONSTRUCTORS AND DESTRUCTOR
  // ---------------------------------------------------------
  /**
   * \brief This constructor creates the object using a Sequence
   * reference and a FullQuality reference.
   *
   * The two references passed to this contructor must remain
   * valid during the entire time of utilization of the instance
   * to be created with the call to this constructor.
   */
  FullyQualifiedSequence(Sequence* seq, FullQuality<C>* qual);

  // ---------------------------------------------------------
  //                    GET AND SET METHODS
  // ---------------------------------------------------------
  /**
   * Return the stored Sequence reference.
   * \return The stored Sequence reference
   */
  Sequence* getSequence() const;
  /**
   * Return the stored Quality reference.
   * \return The stored Quality reference
   */
  FullQuality<C>* getQuality() const;
};


template<class C>
FullyQualifiedSequence<C>::FullyQualifiedSequence(Sequence* seq, FullQuality<C>* qual) {
  this->sequence = seq;
  this->quality = qual;  
}

template<class C>
Sequence* FullyQualifiedSequence<C>::getSequence() const {
  return this->sequence;
}

template<class C>
FullQuality<C>* FullyQualifiedSequence<C>::getQuality() const {
  return this->quality;
}

#endif

