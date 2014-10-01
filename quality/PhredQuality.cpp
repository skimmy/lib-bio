#include "PhredQuality.hpp"

#include <string.h>
#include <math.h>

#include "ProbabilisticQuality.hpp"
#include "../util.h"

#include <iostream>
using namespace std;

#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
using namespace boost;

/*
public:
  
private:

*/


/**************** CONSTRUCTORS AND DESTRUCTOR ***************/

PhredQuality::PhredQuality(const int* quals, size_t n) {
  init(quals,n);
}

PhredQuality::PhredQuality(const string& quals, size_t n) {
  this->n = n;
  this->qualVector = new int[this->n];
  parseQualityString(quals);
}

PhredQuality::PhredQuality(const PhredQuality& other) {
  init(other.qualVector, other.n);
}

PhredQuality::~PhredQuality() {
  if (this->qualVector) {
    delete[] this->qualVector;
  }
  this->qualVector = 0;
}

/**************** 'QUALITY' OVERRIDE METHODS ****************/

double* PhredQuality::getProbabilities(size_t begin, size_t length) const {
  size_t M = ( length != 0 ) ? length : ( this->n - begin);
  size_t N = ( length != 0 ) ? MIN(length, this->n - begin) : M;
  double* out = new double[M];
  for (size_t i = 0; i < N; ++i) {
    out[i] = ProbabilisticQuality::toProbabilistic(this->qualVector[begin+i]);
  }
  return out;
}

int* PhredQuality::getQualities(size_t begin, size_t length) const {
  size_t M = ( length != 0 ) ? length : ( this->n - begin);
  size_t N = ( length != 0 ) ? MIN(length, this->n - begin) : M;
  int* out = new int[M];
  memcpy(out, &(this->qualVector[begin]), N * sizeof(int));
  return out;
}

double PhredQuality::getOverallProbability() const {
  if( this->n <= 0) {
    return -1.0;
  }
  double* probs = getProbabilities();
  double partial = probs[0];
  for (size_t i = 0; i < this->n; ++i) {
    partial *= probs[i];
  }
  delete[] probs;
  return partial;
}

size_t PhredQuality::length() const {
  return this->n;
}

/********************** STATIC METHODS **********************/

int* PhredQuality::toPhred(const double* probs, size_t n) {
  int* quals = new int[n];
  for (size_t i = 0; i < n; ++i) {
    quals[i] = (int)floor( -10 * log10( probs[i]));
  }
  return quals;
}

int PhredQuality::toPhred(double p) {
  return (int)floor(-10 * log10(p));
}



/********************** UTILITY METHODS *********************/

void PhredQuality::init(const int* v, size_t n) {
  this->n = n;
  this->qualVector = new int[this->n];
  memcpy(this->qualVector, v, this->n * sizeof(int));
}

void  PhredQuality::parseQualityString(const string& s) {
  tokenizer<> tok(s);
  size_t i = 0;
  BOOST_FOREACH(string token, tok) {
    qualVector[i++] = atoi(token.c_str());
  }
}

/************************************************************/
