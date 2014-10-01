#include "QualifiedSequence.hpp"

/**************** CONSTRUCTORS AND DESTRUCTOR ***************/

QualifiedSequence::QualifiedSequence(Sequence* seq, Quality* qual) {
  this->sequence = seq;
  this->quality = qual;
}

QualifiedSequence::~QualifiedSequence() {
  // TODO: References should be deallocated outside
  // to preserve consistency (i.e. they are allocated
  // outside
  /*  if(sequence) {
    delete sequence;
  }
  if(quality) {
    delete quality;
    }*/
}

/******************* GET AND SET METHODS ********************/

Quality* QualifiedSequence::getQuality() const {
  return this->quality;
}

Sequence* QualifiedSequence::getSequence() const {
  return this->sequence;
}

void QualifiedSequence::setQuality(Quality* qual) {
  this->quality = qual;
}

void QualifiedSequence::setSequence(Sequence* seq) {
  this->sequence = seq;
}

/************************************************************/
