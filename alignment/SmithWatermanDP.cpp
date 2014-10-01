#include "SmithWatermanDP.hpp"
#include "../util.h"

/***************** CONSTRUCTOR(S)/DESTRUCTOR ****************/

SmithWatermanDP::SmithWatermanDP(const char* s1, size_t n1, const char* s2, size_t n2) 
  // recall that the SW matrix has n+1 and m+1 rows and columns respectively
  : DynamicProgramming(n1+1,n2+1)
{
  this->x = s1;
  this->y = s2;
  // default gap penalty is 1 and efault sim matrix is 1 on
  // diagonal and 0 elsewere
  this->gapPenalty = 1;
  this->sim = createDefaultSimilarityMatrix(256);
  createMatrix();
}

SmithWatermanDP::~SmithWatermanDP() {
  destroyMatrix();
}

/********************* OVERIDDEN METHODS ********************/

void SmithWatermanDP::initMatrix() {
  for (size_t i = 0; i < n; i++ ) {
    matrix[i][0] = 0;
  }
  for (size_t j = 0; j < m; j++) {
    matrix[0][j] = 0;
  }
}

void SmithWatermanDP::computeEntry(size_t i, size_t j) {
  //  int tmp = MAX(1,3);
  int tmp = MAX( (matrix[i-1][j] - gapPenalty), (matrix[i][j-1] - gapPenalty));
  matrix[i][j] = MAX((tmp), (matrix[i-1][j-1] + sim[(int)x[i-1]][(int)y[j-1]]));
}

void SmithWatermanDP::computeMatrix() {
  for (size_t i = 1; i < n; i++ ) {
    for (size_t j = 1; j < m; j++ ) {
      computeEntry(i,j);
    }
  }
}

void SmithWatermanDP::printMatrix() {
  cout << "\t\t";
  for (size_t j = 0; j <m-1; j++ ) {
    cout << y[j] << '\t';
  }
  cout << '\n';
  for (size_t i = 0; i < n; i++ ) {
    if (i > 0) {
      cout << x[i-1];
    }
    cout << '\t';
    for (size_t j = 0; j < m; j++ ) {
      cout << matrix[i][j] << '\t';
    }
    cout << '\n';
  }
}

/*********************** SW PARAMETERS **********************/

void SmithWatermanDP::setSimilarityMatrix(int** s) {
  this->sim = s;
}

void SmithWatermanDP::setGapPenalty(int d) {
  this->gapPenalty = d;
}


/******************* PRIVATE UTIL METHODS *******************/

void SmithWatermanDP::createMatrix() {
  this->matrix = new int*[n];
  for (size_t i = 0; i < n; i++ ) {
    this->matrix[i] = new int[m];
  }
}

void SmithWatermanDP::destroyMatrix() {
  if (this->matrix) {
    for (size_t i = 0; i < n; i++ ) {
      delete[] matrix[i];
    }
    delete[] matrix;
  }
}

int** SmithWatermanDP::createDefaultSimilarityMatrix(size_t s) {
  int** sim = new int*[s];
  for (size_t i = 0; i < s; i++) {
    sim[i] = new int[s];
  }
  for (size_t i = 0; i < s; i++) {
    for (size_t j = 0; j < s; j++ ) {
      sim[i][j] = (i == j) ? 1 : 0;
    }
  }
  return sim;
}
