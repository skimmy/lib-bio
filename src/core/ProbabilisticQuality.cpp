#include <string.h>

#include "../util.h"
#include "../quality.h"

using namespace std;

/*************** CONSTRUCTORS AND DESTRUCTOR ****************/

ProbabilisticQuality::ProbabilisticQuality(const double* prob, size_t n) {
  init(prob,n);
}

ProbabilisticQuality::ProbabilisticQuality(const ProbabilisticQuality& other) {
  init(other.probVector, other.n);
}

ProbabilisticQuality::~ProbabilisticQuality() {
  if (probVector) {
    delete[] probVector;
  }  
  probVector = NULL;
  if (iidProb != NULL) {
    delete iidProb;
  }
  iidProb = NULL;
}

/**************** 'QUALITY' OVERRIDE METHODS ****************/

double* ProbabilisticQuality::getProbabilities(size_t begin, size_t length) const {
  size_t M = ( length != 0 ) ? length : ( this->n - begin );
  size_t N = ( length != 0 ) ? MIN(length, this->n - begin) : M;
  double* out = new double[M];
  memcpy(out, &probVector[begin], N * sizeof(double));
  return out;
}

int* ProbabilisticQuality::getQualities(size_t begin, size_t length) const {
  size_t M = ( length != 0 ) ? length : ( this->n - begin);
  size_t N = ( length != 0 ) ? MIN(length, this->n - begin) : M;
  int* out = new int[M];
  for (size_t i = 0; i < N; ++i) {
    out[i] = PhredQuality::toPhred(this->probVector[begin+i]);
  }
  return out;
}

double ProbabilisticQuality::getOverallProbability() {
  if (iidProb == NULL) {
    iidProb = new double;
    double partial = 1.0 - this->probVector[0];
    for (size_t i = 1; i < this->n; ++i) {
      partial *= (1.0 - probVector[i]);
    }
    *iidProb = partial;
  }
  return (*iidProb);
}

size_t ProbabilisticQuality::length() const {
  return this->n;
}

/********************** STATIC METHODS **********************/

double* ProbabilisticQuality::toProbabilistic(const int* quals, size_t n) {
  double* probs = new double[n];
  for (size_t i = 0; i < n; ++i) {
    probs[i] = qual::PHRED[quals[i]];
  }
  return probs;
}

double ProbabilisticQuality::toProbabilistic(int qual) {
  return qual::PHRED[qual];
}
/********************* FACTORY METHODS **********************/

ProbabilisticQuality ProbabilisticQuality::fromEncodedQuality(const string& v, qual::QualityEncodingType enc) {
  size_t m = v.size();
  double* p = new double[m];
  switch(enc) {
  case qual::QualityEncodingType::SANGER:
    PhredQuality::fromSangerQualities(v,p);
    break;
  case qual::QualityEncodingType::ILLUMINA:
    PhredQuality::fromIlluminaQualities(v,p);
    break;    
  case qual::QualityEncodingType::SOLEXA:
    PhredQuality::fromSolexaQualities(v,p);
    break;    
  case qual::QualityEncodingType::UNKNOWN:
    break;    
  }
//  ProbabilisticQuality pq(p,m);
  delete[] p;
  return pq;
}


/********************* UTILITY METHODS **********************/

void ProbabilisticQuality::init(const double* v, size_t n) {
  this->n = n;
  this->probVector = new double[n];
  memcpy(this->probVector, v, this->n * sizeof(double));
}

/************************************************************/
