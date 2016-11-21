#ifndef UTIL_H
#define UTIL_H

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

size_t editDistanceEncoded(uint64_t s1, size_t n1, uint64_t s2, size_t n2, size_t** dpMatrix);
size_t editDistanceLinSpace(const std::string& s1, const std::string& s2, size_t* v0, size_t* v1);

void printString(char* s, size_t n);
void printDoubleMatrix(double** M, size_t n, size_t m);
double** initDoubleMatrix(size_t n, size_t m);
void clearDoubleMatrix(double** matrix, size_t n, size_t m);
double elementsSumDoubleMatrix(double** matrix, size_t n, size_t m);

Read randomRead(size_t m);

#endif
