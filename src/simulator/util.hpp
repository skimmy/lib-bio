#ifndef UTIL_H
#define UTIL_H

#include <iostream>

#define MAX(X,Y) ( (X) > (Y)  ? (X) : (Y) )
#define MIN(X,Y) ( (X) < (Y)  ? (X) : (Y) )

void initUtil();
void clearUtil();

uint64_t string2Encode(const std::string&s);
std::string encoding2String(uint64_t e, size_t n);

size_t hammingDistance(const char* s1, const char* s2, size_t m);
size_t hammingDistance(const std::string& s1, const std::string& s2, size_t m);
size_t prefixSuffixHammingDistance(const std::string& s1, const std::string& s2, size_t k);
size_t bestHammingOverlap(const std::string& s1, const std::string& s2);

void printString(char* s, size_t n);
void printDoubleMatrix(double** M, size_t n, size_t m);
double** initDoubleMatrix(size_t n, size_t m);
void clearDoubleMatrix(double** matrix, size_t n, size_t m);
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


void print_warning(const std::string &msg);
void fatal_error(const std::string &msg, int exit_code);

#endif
