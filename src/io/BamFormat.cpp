#include "BamFormat.hpp"

using namespace lbiobam;

BamFormat::BamFormat()
  : Format("BAM"), inFilePath("")
{
  
}

void
BamFormat::open(const std::string& filePath)
{
  inFilePath = filePath;
  
  #ifdef HAVE_HTSLIB
  inFile = sam_open(inFilePath.c_str(), "r");
  head = sam_hdr_read(inFile);
  content = bam_init1();
  
  #else
  std::cout << "Error undefined htslib" << std::endl;
  exit(1);
  #endif
}

void
BamFormat::close()
{
  inFilePath = "";
  #ifdef HAVE_HTSLIB  
  if (head)
    {
        bam_hdr_destroy(head);
    }
  if (content)
    {
        bam_destroy1(content);
    }
  if (inFile)
    {
      sam_close(inFile);
      inFile = NULL;
    }
  #endif
}

std::unique_ptr<UIntList>
BamFormat::getAlignmentPositions()
{
  
  std::unique_ptr<UIntList> pList(new UIntList());

  #ifdef HAVE_HTSLIB
  

  while(sam_read1(inFile, head, content) >= 0) {
    pList->push_back(content->core.pos);
  }

  #endif
  
  return pList;
}

IdPos
BamFormat::getNext() {
  IdPos idPos;
  #ifdef HAVE_HTSLIB
  if (sam_read1(inFile, head, content) >= 0) { 
    idPos.first = std::string(bam_get_qname(content));
    idPos.second = content->core.pos;
  }
  else
    {
      idPos.first = "EOF";
      idPos.second = -2;
    }
  #endif
  return idPos;
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
