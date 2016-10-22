#include "BamFormat.hpp"

using namespace lbiobam;

BamFormat::BamFormat()
  : Format("BAM")
{
  
}

void
BamFormat::open(const std::string& filePath)
{
  #ifdef HAVE_HTSLIB
  inFile = sam_open(filePath.c_str(), "r");
  #elif
  std::cout << "Error undefined htslib" << std::endl;
  exit(1);
  #endif
}

void
BamFormat::close()
{
  #ifdef HAVE_HTSLIB
  if (inFile)
    {
      sam_close(inFile);
      inFile = NULL;
    }
  #elif
  #endif
}

std::unique_ptr<UIntList>
BamFormat::getAlignmentPositions()
{
  std::unique_ptr<UIntList> pList(new UIntList());
  bam_hdr_t* head = sam_hdr_read(inFile);
  bam1_t* content = bam_init1();
  while(sam_read1(inFile, head, content) >= 0) {
    pList->push_back(content->core.pos);
  }
  bam_destroy1(content);
  bam_hdr_destroy(head);
  return pList;
}

std::string
BamFormat::loadFromFile(const std::string& fileName)
{
  open(fileName);
  return "";
}

std::string
BamFormat::getSequence() const
{
  return "";
}

std::string
BamFormat::getHeader() const
{
  return "";
}
