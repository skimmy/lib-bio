#include "FastqRead.hpp"

#include <core/ProbabilisticQuality.hpp>

#include <fstream>

// Remove when 'getFileLength' will be moved to new util sub-project
#include <sys/stat.h>

size_t getFileLength(const std::string& filePath) {
  size_t n = 0;
  struct stat buffer;
  int status = stat(filePath.c_str(), &buffer);
  if (status == 0) {
    n = (size_t) buffer.st_size;
  }
  return n;
}


#define MAX_BUFFER_SIZE 4096

/******************** SUPPORT FUNCTIONS *********************/
size_t nextReadPosition(std::ifstream& ifs, size_t off = 0) {
  ifs.seekg(off);
  size_t candidate = off;
  char line[MAX_BUFFER_SIZE];
  line[0] = 0;
  // look for the first occurance of '@' at the beginning of a line
  while(line[0] != '@' && !ifs.eof()) {
    candidate = ifs.tellg();
    ifs.getline(line, MAX_BUFFER_SIZE);
  }
  // end of file reached and no reads have been found
  if (ifs.eof()) {
    return off;
  }
  size_t tmp = ifs.tellg();  
  // read the next line to double check that '@' defines a proper start of a read
  ifs.getline(line, MAX_BUFFER_SIZE);
  candidate = (line[0] != '@') ? candidate : tmp;
  return candidate;
}


/********************** CONSTRUCTOR(S) **********************/

FastqRead::FastqRead() 
  : Read()
{
  this->quality = NULL;
}

FastqRead::FastqRead(const FastqRead& other) {
  bases = other.bases;
  qualities =other.qualities;
  header = other.header;
  this->quality = NULL;
}

FastqRead::~FastqRead() {
  delete this->quality;
  this->quality = NULL;
}

//////////////////////////////////////////////////////////////////////
//                      STREAM OPERATORS
//////////////////////////////////////////////////////////////////////

std::istream& operator>>(std::istream& is, FastqRead& read) {

  // This works as long as the fastq is well formed. More robust code
  // should be implemented in the future. For example aborting current
  // read when malformed file is encountered.

  std::getline(is, read.header);
  std::getline(is, read.bases);
  // The '+' divider between sequence and qualities is ignored
  std::string trash_string;
  std::getline(is, trash_string);
  std::getline(is, read.qualities);

  return is;
}

std::ostream& operator<<(std::ostream& os, const FastqRead& read) {
  os << read.header << '\n';
  os << read.bases << '\n';
  os << "+\n";
  os << read.qualities << '\n';
  return os;
}

/********************** STATIC METHODS **********************/

std::vector< size_t > FastqRead::splitReads(const std::string& filePath, size_t T) {
  std::vector< size_t > offsets(T);  
  std::ifstream ifs(filePath);
  size_t n = getFileLength(filePath);
  size_t step = n / T;
  offsets[0] = 0;
  size_t i = 1;
  while(i < T) {
    offsets[i] = nextReadPosition(ifs, offsets[i-1] + step);
    ++i;
  }
  return offsets;
}

/************************************************************/

void FastqRead::autoDecodeQualities() {
  this->quality = NULL;
  // if no quality string is available no Quality object is stored (i.e., NULL)
  if (!hasQualities()) {    
    return;
  }
  // try to 'guess' the encoding used
  qual::QualityEncodingType encoding = qual::parseQuality(this->qualities);
  this->quality = new ProbabilisticQuality(ProbabilisticQuality::fromEncodedQuality(this->qualities, encoding));
  
}
