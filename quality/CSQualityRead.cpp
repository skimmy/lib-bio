#include "CSQualityRead.hpp"

#include <iostream>
#include <boost/tokenizer.hpp>
using namespace std;
using namespace boost;

/**************** CONSTRUCTORS AND DESTRUCTOR ***************/

CSQualityRead::CSQualityRead(size_t n, const string& qualities) 
  : ReadQuality(n) 
{
  computeValues(qualities);
}

CSQualityRead::CSQualityRead(const CSQualityRead& other)
  : ReadQuality(other)
{
}

/***************** PROBABILITIES COMPUTATION ****************/

void CSQualityRead::computeValues(const string& quals) {
  tokenizer<> tok(quals);
  size_t i = 0;
  for (tokenizer<>::iterator beg=tok.begin(); beg!=tok.end(); beg++) {
    this->qualities[i] = atoi((*beg).c_str());
    i++;
  }
  computeProbabilities();
  for (size_t j = 0; j < length; ++j) {
    cout << probabilities[j] << "(" << (int)qualities[j] << ")" << "\n";
  }
  cout << '\n';
}

void CSQualityRead::computeProbabilities() {
  for (size_t i = 0; i < length; ++i) {
    probabilities[i] = probLookup[qualities[i]];
  }
}

/************************************************************/
