#include "MatchSimilarity.hpp"

/**************** CONSTRUCTOR(S) DESTRUCTOR ******************/

MatchSimilarity::MatchSimilarity(string* s1, string* s2) {
  this->seq1 = s1;
  this->seq2 = s2;
}

MatchSimilarity::~MatchSimilarity() {
  if (this->matrix) {
    for (size_t i = 0; i < n; ++i) {
      delete[] this->matrix[i];
    }
    delete[] this->matrix;
    this->matrix = 0;
  }
}

/************** DYNAMICPROGRAMMING OVERIDE ******************/

void MatchSimilarity::initMatrix() {
}

void MatchSimilarity::computeEntry(size_t i, size_t j) {
}

void MatchSimilarity::printMatrix() {
}

/************************************************************/
