#ifndef UTIL_H
#define UTIL_H

#include <iostream>

// general variables
// extern char bases[];
// extern char revBases[128];
extern double* power4_lookup;


void initUtil(size_t m);
void clearUtil();

void fatal_error(const std::string &msg, int exit_code);

uint64_t string2Encode(const std::string&s);
std::string encoding2String(uint64_t e, size_t n);

size_t hammingDistance(const char* s1, const char* s2, size_t m);
size_t hammingDistance(const std::string& s1, const std::string& s2, size_t m);
size_t prefixSuffixHammingDistance(const std::string& s1, const std::string& s2, size_t k);
size_t bestHammingOverlap(const std::string& s1, const std::string& s2);

//////////////////////////////////////////////////////////////////////
//
//         FUNCTIONS FOR OPERATING ON VECTOR AND MATRIXES
//
//////////////////////////////////////////////////////////////////////

// TODO: can be templated
double elementsSumDoubleMatrix(double** matrix, size_t n, size_t m);


template<typename T>
void printVector(T* v, size_t n, char* delim = " ") {
  for (int i =0; i <n;i++){
    std::cout << v[i] << delim;
  }
  std::cout << std::endl;
}

template<typename T>
void printMatrix(size_t n, size_t m, T** mat, std::string delim = " ") {
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      std::cout << mat[i][j] << delim;
    }
    std::cout << "\n";
  }
}

template<typename T>
T* allocVector(size_t n) {
  return new T[n];
}

template<typename T>
void freeVector(T* v) {
  delete[] v;
}

template<typename T>
T** allocMatrix(size_t n, size_t m) {
  T** mat = new T*[n];
  for (size_t i = 0; i < n; ++i) {
    mat[i] = allocVector< T >(m);
  }
  return mat;
}

template<typename T>
void freeMatrix(size_t n, size_t m, T** mat) {
  for (size_t i = 0; i < n; ++i) {
    delete[] mat[i];
  }
  delete[] mat;
}

template<typename T>
void resetMatrix(size_t n, size_t m, T** mat, T val) {
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      mat[i][j] = val;
    }
  }
    
}


template<typename T>
void writeVectorOnStream(T v[], size_t n, std::ostream& os, std::string delim = "\n") {
  for (size_t i = 0; i < n; ++i) {
    os << v[i] << delim;
  }
  os.flush();
}

//////////////////////////////////////////////////////////////////////
//
//                      MATH RELATED FUNCTIONS
//
//////////////////////////////////////////////////////////////////////

// Here there are some useful functions and clases for math operations
// and property testing. So far they are few, if their number
// increases consider moving to a dedicated header file.


template<typename T>
class GeometricProgression {
public:
  GeometricProgression(double r, T s)
    : ratio(r), start(s), next(s)
  {  }
  
  T getNext() {
    T tmp =next;
    next *= ratio;
    return tmp; 
  }

  T getCurrent() {
    return next;
  }

  void rewind() {
    next = start;
  }

  std::vector<T> valuesLeq(T maxVal) {
    std::vector<T> vals;
    rewind();
    while(next <= maxVal) {
      vals.push_back(getNext());
    }
    return vals;
  }
  
private:
  const double ratio;
  const T start;
  T next;
};

template<typename T>
class LinearProgression {
public:
  LinearProgression(T s, T b)
    : step(s), start(b), next(b) {}
  
  T getNext() {
    T tmp = next;
    next += step;
    return tmp;
  }

  T getCurrent() {
    return next;
  }

  void rewind() {
    next = start;
  }

  std::vector<T> valuesLeq(T maxVal) {
    std::vector<T> vals;
    rewind();
    while(next <= maxVal) {
      vals.push_back(getNext());
    }
    return vals;
  }
  
private:
  const T step;
  const T start;
  T next;

};


// monotonicity tests
template<typename T>
bool isMonotoneNonDecreasing(const T& v, lbio_size_t n) {
  for (lbio_size_t i = 1; i < n; ++i) {
    if (v[i] < v[i-1]) {
      return false;
    }
  }
  return true;
}

template<typename T>
bool isMonotoneNonIncreasing(const T& v, lbio_size_t n) {
  for (lbio_size_t i = 1; i < n; ++i) {
    if (v[i] > v[i-1]) {
      return false;
    }
  }
  return true;
}

template<typename T>
bool isConstant(const T& v, lbio_size_t n) {
  for (lbio_size_t i = 0; i < n; ++i) {
    if (v[i] != v[0]) {
      return false;
    }
  }
  return true;  
}

// elementwise equality test
template<typename T>
bool areElementwiseEqual(const T& v1, const T& v2, lbio_size_t n) {
  for (lbio_size_t i = 0; i < n; ++i) {
    if (v1[i] != v2[i]) {
      return false;
    }
  }
  return true;
}

#endif
