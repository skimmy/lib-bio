#ifndef POSITION_H
#define POSITION_H

#include <iostream>

/**
 * \brief This class represents a \e position within a set of
 * sequences.
 *
 * The position is in this class threated with a wider meaning,
 * given a set of sequences a <em>position in the</em> is given
 * by the \e id of the sequence within the set and by the position
 * in the sequnce.
 * While the id could be (in principle) an arbitrary \e comparable
 * object (i.e. implementing the comparison operators <tt>==, !=)</tt>
 * and assignable (i.e. they must implement <tt>==</tt>
 * operator and a copy constructor) the position is of \c size_t
 * type to reflect the fact that we really are working with position
 */
template<class ID> 
class Position {
protected:
  ID id;
  size_t position;
public:
  // ---------------------------------------------------------
  //                 CONSTRUCTORS AND DESTRUCTOR
  // ---------------------------------------------------------
  /**
   * \brief Creates an empty Position object.
   *
   * The content of the object should not be considered as long
   * as proper set methods has been called because (for efficiency
   * and for templating reasons) they are noti initialized
   * \sa Position(ID i, size_t p)
   * \sa Position(const Position& other)
   */
  Position();
  /**
   * \brief Creates a Position object with the given content
   *
   * \param i The \e id object
   * \param p The position
   *
   * \sa Position()
   * \sa Position(const Position<ID>& other)
   */
  Position(ID i, size_t p);
  /**
   * \brief Creates a copy of the given Position object
   *
   * \param other The Position object to be copied
   *
   * \sa Position(ID i, size_t p)
   * \sa Position()
   * \sa operator=()
   */
  Position(const Position<ID>& other);
  
  // ---------------------------------------------------------
  //                    GET AND SET METHODS
  // ---------------------------------------------------------
  /**
   * \brief Returns the sequence \e id
   *
   * \return The sequence \e id
   *
   * \sa setSequenceId()
   * \sa getPosition()
   */
  ID getSequenceId() const;
  /**
   * \brief Returns the \e position
   *
   * \return The \e position
   *
   * \sa setPosition()
   * \sa getSequenceId()
   */
  size_t getPosition() const;
  /**
   * \brief Sets the \e id
   * 
   * \param i The new \e id
   *
   * \sa getSequenceId()
   * \sa setPosition()
   */
  void setSequenceId(ID i);
  /**
   * \brief Sets the \e position
   *
   * \param p The new \e position
   *
   * \sa getPosition()
   * \sa setSequenceId()
   */
  void setPosition(size_t p);
  
  // ---------------------------------------------------------
  //                        OPERATORS
  // ---------------------------------------------------------
  /**
   * \brief Tests equality between two Position objects
   *
   * The test is performed on both components (i.e. \e id
   * and \e position)
   *
   * \return \c true if the two objects are equal, \c false otherwise
   *
   * \sa operator!=()
   */
  bool operator==(const Position& other);
  /**
   * \brief Tests inequality between two Position objects
   *
   * The test is performed on both components (i.e. \e id 
   * and \e position)
   *
   * \return \c true if the two objects are different, \c false otherwise
   *
   * \sa operator==()
   */
  bool operator!=(const Position& other);

  /**
   * \brief Assign a Poisiton object to another
   *
   * \param other The Position to be assigned
   * \return The object after assignment
   *
   * \sa Position(const Position<ID>& other)
   */
  Position& operator=(const Position& other) {
    this->id = other.id;
    this->position = other.position;
    return *this;
  } 

};

template<class ID> 
Position<ID>::Position() {  }

template<class ID> 
Position<ID>::Position(ID i, size_t p) {
  this->id = i;
  this->position = p;
}

template<class ID> 
Position<ID>::Position(const Position<ID>& other) {
  this->id = other.id;
  this->position = other.position;
}

template<class ID> 
ID Position<ID>::getSequenceId() const {
  return this->id;
}

template<class ID>
size_t Position<ID>::getPosition() const {
  return this->position;
}

template<class ID> 
void Position<ID>::setSequenceId(ID i) {
  this->id = i;
}

template<class ID> 
void Position<ID>::setPosition(size_t p) {
  this->position = p;
}

template<class ID> 
bool Position<ID>::operator==(const Position& other) {
  return ( (this->id == other.id) && (this->position == other.position) );
}

template<class ID> 
bool Position<ID>::operator!=(const Position& other) {
  return !( (this->id == other.id) && (this->position == other.position) );
}



#endif
