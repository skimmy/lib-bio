#ifndef KEY_VALUE_PAIR_H
#define KEY_VALUE_PAIR_H

/**
 * \brief This class represents a \e pair of two elements.
 *
 * The first element is used as \e key while the second element
 * is used as \e value. This distinction is important for those
 * classes wanting to be used as template parameters. In particular
 * keys must implement the <em>equality operator</em> <tt>= =</tt>
 *
 * For efficiency reason KeyValuePair uses pointers rather
 * than 'real' objects, this means that all the reference
 * passed to the object must point to a valid instance of
 * templated classes during the entire scope of the instance
 *
 */
template<class K, class V>
class KeyValuePair {
private:
  K* key;
  V* value;
public:
  // ---------------------------------------------------------
  //                 CONSTRUCTORS AND DESTRUCTOR
  // ---------------------------------------------------------
  /**
   * \brief Constructs a KeyValuePair object with the given
   * key and value elements
   *
   * \param k The key element of the pair
   * \param v The value element of the pair
   */
  KeyValuePair(K* k, V* v);
  // ---------------------------------------------------------
  //                    GET AND SET METHODS
  // ---------------------------------------------------------
  /**
   * \brief Returns the key element of the pair
   *
   * \return The key element
   *
   * \sa setKey()
   * \sa getValue()
   */
  K* getKey() const;
  /**
   * \brief Returns the value element of the pair
   *
   * \return The value element
   *
   * \sa setValue()
   * \sa getKey()
   */
  V* getValue() const;
  /**
   * \brief Sets the key element of the pair
   *
   * \param k The new key
   *
   * \sa getKey()
   * \sa setValue()
   */
  void setKey(K* k);
  /**
   * \brief Sets the value element of the pair
   *
   * \param v The new value
   *
   * \sa getValue()
   * \sa setKey()
   */
  void setValue(V* v);
  // ---------------------------------------------------------
  //                         OPERATORS
  // ---------------------------------------------------------
  /**
   * \brief Tests the equality between two KeyValuePair objects
   *
   * The test is performed only on the key element, this is
   * reasonable bearing in mind that key value pairs are usually
   * used to index an object (i.e. the value) using another
   * object (i.e. the key) as index.
   *
   * \param other The KeyValuePair to compare with the actual object
   * \return \c true if the two KeyValuePair are equal, \c false otherwise
   *
   * \sa operator!=()
   * 
   */
  bool operator == (const KeyValuePair<K,V>& other) const;
  /**
   * \brief Test the inequality between two KeyValuePair objects
   *
   * The test is performed only on the key element, this is
   * reasonable bearing in mind that key value pairs are usually
   * used to index an object (i.e. the value) using another
   * object (i.e. the key) as index.
   * 
   * \param other The KeyValuePair to compare with the actual object
   * \return \c false if the two KeyValuePair are equal, \c true otherwise
   *
   * \sa operator==()
   * 
   */
  bool operator != (const KeyValuePair<K,V>& other) const;
};

template< class K, class V >
KeyValuePair<K,V>::KeyValuePair(K* k, V* v) {
  this->key = k;
  this->value = v;
}

template< class K, class V >
K* KeyValuePair<K,V>::getKey() const {
  return this->key;
}

template< class K, class V >
V* KeyValuePair<K,V>::getValue() const {
  return this->value;
}

template< class K, class V >
void KeyValuePair<K,V>::setKey(K* k) {
  this->key = k;
}

template< class K, class V >
void KeyValuePair<K,V>::setValue(V* v) {
  this->value = v;
}

template< class K, class V >
bool KeyValuePair<K,V>::operator == (const KeyValuePair<K,V>& other) const {
  return ( this->key == other.key );
}

template< class K, class V >
bool KeyValuePair<K,V>::operator != (const KeyValuePair<K,V>& other) const {
  return !( this->key == other.key );
}

#endif

