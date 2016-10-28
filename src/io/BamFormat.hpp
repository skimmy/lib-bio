#ifndef BAMFORMAT_H
#define BAMFORMAT_H

#include "../io.h"

#include <list>
#include <memory>

typedef std::list<uint64_t> UIntList;
typedef std::pair<std::string,int64_t> IdPos;

// This is a variable to be set using autoconf & co
#define HAVE_HTSLIB 1

#ifdef HAVE_HTSLIB

#include <htslib/sam.h>

#endif

namespace lbiobam
{
  /**
   * This <code>struct</code> contains the fields than can be found
   * in a BAM files. 
   *
   * This struct is inspired by
   * <code>bam1_core_t</code> and <code>bam1_t</code> from the
   * <code>htslib</code> package of <b>Samtools</b>. The original
   * structures from the library are not used here to maintain
   * compatibility in the case of missing library.
   */
  typedef struct
  {
    
    int32_t tid;
    int32_t pos;
    uint32_t bin:16, qual:8, l_qname:8;
    uint32_t flag:16, n_cigar:16;
    int32_t l_qseq;
    int32_t mtid;
    int32_t mpos;
    int32_t isize;

    int l_data, m_data;
    uint8_t *data;
     
    
  } BamAlignInfo;
  
  enum BamOpenMode { BamOpenRead, BamOpenWrite, BamOpenAppend, NotOpened };

  class BamFormat : public Format {
  private:
    #ifdef HAVE_HTSLIB
    htsFile* hFile = NULL;
    bam_hdr_t* head = NULL;
    bam1_t* content = NULL;

    #endif

    BamOpenMode mode;
    std::string hFilePath;
  public:
    BamFormat();
    void open(const std::string& filePath, BamOpenMode opMode = BamOpenRead);
    void close();

    std::unique_ptr<UIntList> getAlignmentPositions();
    IdPos getNext();

    void setBamHeader();
    void copyHeader(const BamFormat& other);
    void writeBamHeader();

    // As of now these are meaningless and should not be used.
    // In the future we may give them speial meaning (bad design [sic]).
    std::string loadFromFile(const std::string& fileName) { return ""; }
    std::string getSequence() const  { return ""; }
    std::string getHeader() const  { return ""; }
  };
  
}

#endif
