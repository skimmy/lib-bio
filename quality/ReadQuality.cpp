#include "ReadQuality.hpp"

#include <math.h>
#include <string.h>

double probLookup[256];

/**************** CONSTRUCTORS AND DESTRUCTOR ***************/

ReadQuality::ReadQuality(size_t n) {
  init(n);
}

ReadQuality::ReadQuality(const ReadQuality& other) {
  this->length = other.length;
  this->qualities = new uint8_t[this->length];
  this->probabilities = new double[this->length];
  memcpy(this->qualities, other.qualities, this->length * sizeof(uint8_t));
  memcpy(this->probabilities, other.probabilities, this->length * sizeof(double));
}

ReadQuality::~ReadQuality() {
  if(this->qualities) {
    delete[] this->qualities;
    this->qualities = 0;
  }

  if(this->probabilities) {
    delete[] this->probabilities;
    this->probabilities = 0;
  }

  this->length = -1;
}

/****************** PRIVATE UTILITY METHODS *****************/

void ReadQuality::init(size_t n) {
  this->length = n;
  this->qualities = new uint8_t[n];
  this->probabilities = new double[n];
}

/******************* GET AND SET METHODS ********************/

size_t ReadQuality::getLength() const {
  return this->length;
}

uint8_t* ReadQuality::getQualities() {
  return this->qualities;
}

double* ReadQuality::getProbabilities() {
  return this->probabilities;
}

double ReadQuality::getOverallProbability() const {
  if (!probabilities || length <= 0) {
    return -1.0;
  }
  double x = probabilities[0];
  for (size_t i = 0; i < length; ++i) {
    x *= probabilities[i];
  }
  return x;
}

/************************************************************/

void initLookupTables() {
  for (size_t i = 0; i < 256; ++i) {
    probLookup[i] = pow(10.0,  -1 * (double)(i) / 10.0);
  }
}




