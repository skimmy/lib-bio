#ifndef _HYBRID_INDEX_
#define _HYBRID_INDEX_

template<class T, class C>
class HybridIndex {
private:
  C* index;
public:
  void add(T& item);
};

#endif

template<class T, class C>
void HybridIndex<T,C>::add(T& item) {
}
