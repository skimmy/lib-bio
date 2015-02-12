#ifndef _HYBRID_INDEX_
#define _HYBRID_INDEX_

#include <unordered_map>
#include <cstdlib>

#include <iostream>

template<class K, class V>
class HybridIndex {
private:
  std::unordered_map<K, V>* index;
  size_t k_prime;
public:
  HybridIndex();
  ~HybridIndex();
  void add(K& key, V& value);
};

#endif

template<class K, class V>
HybridIndex<K,V>::HybridIndex() {
  k_prime = 1 << 8; // TODO: Change to be dynamically tuned
  index = new std::unordered_map<K, V>[k_prime];
}

template<class K, class V>
HybridIndex<K,V>::~HybridIndex() {
  if (index != NULL) {
    delete[] index;
  }
}
template<class K, class V>
void HybridIndex<K,V>::add(K& key, V& value) {
  size_t i = (key % k_prime);
  this->index[i][key] = value;
  std::cout << this->index[i].size() << std::endl;
}

