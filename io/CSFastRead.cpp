#include "CSFastRead.hpp"

#define MAX_BUFFER_SIZE 4096

// ---------------------------------------------------------
//                      CONSTRUCTORS
// ---------------------------------------------------------

CSFastRead::CSFastRead()
  : Read(), qualityFileName("")
{
   
}
 
CSFastRead::CSFastRead(const CSFastRead& other) {
  bases = other.bases;
  qualities = other.qualities;
  header = other.header;
  primer = other.primer;
}

// ---------------------------------------------------------
//                    GET AND SET METHODS
// ---------------------------------------------------------

char CSFastRead::getPrimer() const {
  return this->primer;
}

void CSFastRead::setPrimer(char p) {
  this->primer = p;
}

list< FullyQualifiedSequence< ColorAlphabet > > CSFastRead::getFullyQualifiedKMerList(size_t k) const {
  list< FullyQualifiedSequence< ColorAlphabet > > qualKmersList;
  size_t N = bases.size() - k + 1;
  // compute raw string and quality for entire read
  const char* rawSequence = this->bases.c_str();
  PhredQuality quals(this->qualities, this->bases.size());
  for (size_t i = 0; i < N; ++i) {
    // get quality and k-mer for subsequence [i .. i+k-1] of the read
    PhredQuality qual(quals.getQualities(i,k), k);
    KMer* kmer = new KMer(&rawSequence[i], k);
    // creates and insert into the output list a new entry
    FullQuality< ColorAlphabet >* quality = new FullQuality< ColorAlphabet >(qual,*kmer);
    qualKmersList.push_back(FullyQualifiedSequence< ColorAlphabet >(kmer, quality));
  }
  return qualKmersList;
}

// ---------------------------------------------------------
//                  '>>' AND '<<' OPERATORS
// ---------------------------------------------------------

istream& operator>>(istream& is, CSFastRead& read) {
  // WARNING: this methods doesn't load quality values (yet)
  char buffer[MAX_BUFFER_SIZE];
  // header starts with '>'
  is.getline(buffer,MAX_BUFFER_SIZE);
  read.header = string(buffer);
  // read actual colors
  is.get(read.primer);
  is.getline(buffer,MAX_BUFFER_SIZE);
  read.bases = string(buffer);  
  return is;
}

ostream& operator<<(ostream& os, CSFastRead& read) {
  // WARNING: this methods doesn't write quality values (yet)
  os << read.header << '\n';
  os << read.primer << read.bases << '\n';
  return os;
}

void CSFastRead::loadBasesAndQualitiesFromFiles(istream& bases_is, istream& qual_is) {
  bases_is >> (*this);
  char buffer[MAX_BUFFER_SIZE];
  // discard header
  qual_is.getline(buffer,MAX_BUFFER_SIZE);
  if (string(buffer) != this->header) {
    cerr << "[ERROR] - Inconsistency between colors and qualities header\n";
  }
  // read qualitis string
  qual_is.getline(buffer,MAX_BUFFER_SIZE);
  this->qualities = string(buffer);
}
