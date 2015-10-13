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
  // defualt value for backtracking is not enbaled
  this->btEnabled = false;
  this->btMatrix = NULL;
  createMatrix();
}

SmithWatermanDP::SmithWatermanDP(const string& s1, const string& s2)
  : DynamicProgramming(s1.size() + 1, s2.size() + 1){
  this->x = s1.c_str();
  this->y = s2.c_str();
  // default gap penalty is 1 and efault sim matrix is 1 on
  // diagonal and 0 elsewere
  this->gapPenalty = 1;
  this->sim = createDefaultSimilarityMatrix(256);
  // defualt value for backtracking is not enbaled
  this->btEnabled = false;
  this->btMatrix = NULL;
  createMatrix();
}

SmithWatermanDP::~SmithWatermanDP() {
  destroyMatrix();
  deleteBtMatrix();
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
  // if enabled init the backtrack matrix
  if (this->btEnabled) {
    // init or reset backtrack matrix
    if (this->btMatrix == NULL) {
      this->initBtMatrix();
    } else {
      this->resetBtMatrix();
    }
  }
  // compute the score matrix
  for (size_t i = 1; i < n; i++ ) {
    for (size_t j = 1; j < m; j++ ) {
      computeEntry(i,j);
      // IMPORTANT: this operation must be performed AFTER computeEntry()
      if (this->btEnabled) {
	this->btMatrix[i][j] = getBtForPosition(i,j);
      }
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

/**************** RESEULT RETRIEVING METHODS ****************/
MatrixPoint2D SmithWatermanDP::getGlobalBest() const {
  MatrixPoint2D entryIndices;
  this->indicesOfMaxElement(entryIndices.i, entryIndices.j);
  return entryIndices;
}

int SmithWatermanDP::getScoreAt(const MatrixPoint2D& p) {
  return this->matrix[p.i][p.j];
}


/******************* BACKTRACKING METHODS *******************/

void SmithWatermanDP::enableBacktrack() {
  this->btEnabled = true;
}

void SmithWatermanDP::disableBacktrack() {
  this->btEnabled = false;
}

bool SmithWatermanDP::isBacktrackEnabled() const {
  return this->btEnabled;
}

void SmithWatermanDP::printBacktrackMatrix() const {
  if (!this->btEnabled) {
    cout << "Backtrack not enabled!" << endl;
    return;
  }
  string symbols = "UMSID";
  for (size_t i = 0; i < n; i++ ) {
    cout << " \t";
    for (size_t j = 0; j < m; j++ ) {
      cout << symbols[(int)btMatrix[i][j]] << '\t';
    }
    cout << '\n';
  }
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
    for (size_t i = 0; i < n; ++i) {
      delete[] matrix[i];
    }
    delete[] matrix;
  }
}

int** SmithWatermanDP::createDefaultSimilarityMatrix(size_t s) {
  int** sim = new int*[s];
  for (size_t i = 0; i < s; ++i) {
    sim[i] = new int[s];
  }
  for (size_t i = 0; i < s; ++i) {
    for (size_t j = 0; j < s; ++j) {
      sim[i][j] = (i == j) ? 1 : 0;
    }
  }
  return sim;
}

void SmithWatermanDP::indicesOfMaxElement(size_t& i, size_t& j) const {
  // init with first element of the matrix
  i = j = 0;
  for (size_t r = 0; r < this->n; ++r) {
    for (size_t c = 0; c < this->m; ++c) {
      if (this->matrix[r][c] > this->matrix[i][j]) {
	i = r;
	j = c;
      }
    }
  }
}

void SmithWatermanDP::initBtMatrix() {
  this->btMatrix = new BacktrackOperation*[n];
  for (size_t i = 0; i < n; ++i) {
    this->btMatrix[i] = new BacktrackOperation[m];
    for (size_t j = 0; j < m; ++j) {
      // initially every entry is 'Unset'
      this->btMatrix[i][j] = Unset;      
    }
  }
}

void SmithWatermanDP::deleteBtMatrix() {
  if (this->btMatrix != NULL) {
    for (size_t i = 0; i < n; ++i) {
      delete[] this->btMatrix[i];
    }
    delete[] this->btMatrix;
    this->btMatrix = NULL;
  }  
}

void SmithWatermanDP::resetBtMatrix() {
  for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < m; ++j) {
        this->btMatrix[i][j] = Unset;      
      }
  }
}

BacktrackOperation SmithWatermanDP::getBtForPosition(size_t i, size_t j) const {
  /*
  int tmp = MAX( (matrix[i-1][j] - gapPenalty), (matrix[i][j-1] - gapPenalty));
  matrix[i][j] = MAX((tmp), (matrix[i-1][j-1] + sim[(int)x[i-1]][(int)y[j-1]]));
   */
  // match or substitution
  if (matrix[i][j] == (matrix[i-1][j-1] + sim[(int)x[i-1]][(int)y[j-1]])) {
    if (x[i-1] == y[i-1]) {
      return Match;
    } else {
      return Substitution;
    }
  }
  // deletion
  if (matrix[i][j] == matrix[i][j-i]){
    return Deletion;
  }
  // insertion
  return Insertion;
}
