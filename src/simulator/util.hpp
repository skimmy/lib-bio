#ifndef UTIL_H
#define UTIL_H

void initUtil();
void clearUtil();

size_t hammingDistance(const char* s1, const char* s2, size_t m);
size_t hammingDistance(const std::string& s1, const std::string& s2, size_t m);
size_t prefixSuffixHammingDistance(const std::string& s1, const std::string& s2, size_t k);
size_t bestHammingOverlap(const std::string& s1, const std::string& s2);

void printString(char* s, size_t n);
void printDoubleMatrix(double** M, size_t n, size_t m);
double** initDoubleMatrix(size_t n, size_t m);
void clearDoubleMatrix(double** matrix, size_t n, size_t m);
double elementsSumDoubleMatrix(double** matrix, size_t n, size_t m);

Read randomRead(size_t m);

#endif
