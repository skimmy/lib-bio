#ifndef BAMFORMAT_H
#define BAMFORMAT_H

#include "../io.h"

#include <list>
#include <memory>

typedef std::list<uint64_t> UIntList;

// This is a variable to be set using autoconf & co
#define HAVE_HTSLIB 1

#include <htslib/sam.h>

namespace lbiobam
{

  class BamFormat : public Format {
  private:
    #ifdef HAVE_HTSLIB
    htsFile* inFile = NULL;
    #endif
  public:
    BamFormat();
    void open(const std::string& filePath);
    void close();

    std::unique_ptr<UIntList> getAlignmentPositions();
    
    std::string loadFromFile(const std::string& fileName);
    std::string getSequence() const;
    std::string getHeader() const;
  };
  
}

#endif
