#ifndef SCORED_POSITION_H
#define SCORED_POSITION_H

#include "Position.hpp"

template<class ID, class S>
class ScoredPosition : public Position<ID> {
private:
  S score;
public:  
  ScoredPosition(ID i, size_t p, S s) : Position<ID>(i,p) { this->score = s; }
  S getScore() const { return score; }
};

#endif 
