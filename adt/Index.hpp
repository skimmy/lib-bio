#ifndef INDEX_H
#define INDEX_H

#include "KeyValuePair.hpp"

#include <boost/unordered_map.hpp>
using namespace boost;

/**
 * 
 */
template<class K, class V>
class Index {
private:
  unordered_map<K,V> map;
public:
  void insert(const KeyValuePair<K,V>& pair) { map[*(pair.getKey())] = *(pair.getValue()); }
  // virtual KeyValuePair<K,V>* find(const KeyValuePair& pair); 
  // virtual KeyValuePair<K,V> remove(const KeyValuePair& pair);
};

#endif
