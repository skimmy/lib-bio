#include "BamFormat.hpp"

using namespace lbiobam;

BamFormat::BamFormat()
  : Format("BAM"), hFilePath(""), mode(NotOpened)
{
  
}

void
BamFormat::open(const std::string& filePath, BamOpenMode opMode)
{
  hFilePath = filePath;

  // Read only case
  if (opMode == BamOpenRead)
    {
      mode = BamOpenRead;
      std::cout << hFilePath << " -->  READ\n";
#ifdef HAVE_HTSLIB

      hFile = sam_open(hFilePath.c_str(), "r");
      head = sam_hdr_read(hFile);
      content = bam_init1();
  
#else
      std::cout << "Error undefined htslib" << std::endl;
      exit(1);
#endif
    }

  // (Over)write only case
  if (opMode == BamOpenWrite)
    {
      mode = BamOpenWrite;
      std::cout << hFilePath << " -->  WRITE\n";      
#ifdef HAVE_HTSLIB
      hFile = sam_open(hFilePath.c_str(), "w");
      
#else
      std::cout << "Error undefined htslib" << std::endl;
      exit(1);
#endif
    }
}

void
BamFormat::close()
{
  hFilePath = "";
  #ifdef HAVE_HTSLIB  
  if (head)
    {
        bam_hdr_destroy(head);
    }
  if (content)
    {
        bam_destroy1(content);
    }
  if (hFile)
    {
      sam_close(hFile);
      hFile = NULL;
    }
  #endif
}

std::unique_ptr<UIntList>
BamFormat::getAlignmentPositions()
{
  
  std::unique_ptr<UIntList> pList(new UIntList());
  if (mode != BamOpenRead)
    {     
      return NULL;
    }

  #ifdef HAVE_HTSLIB
  

  while(sam_read1(hFile, head, content) >= 0) {
    pList->push_back(content->core.pos);
  }

  #endif
  
  return pList;
}

IdPos
BamFormat::getNext() {
  IdPos idPos;
  if (mode != BamOpenRead)
    {
      idPos.first = "ERR";
      idPos.second = -10;
      return idPos;
    }
  #ifdef HAVE_HTSLIB
  if (sam_read1(hFile, head, content) >= 0) { 
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


void
BamFormat::setBamHeader()
{  
}

void
BamFormat::copyHeader(const BamFormat& other)
{
#ifdef HAVE_HTSLIB
  
  this->head = bam_hdr_dup(other.head);

#else
    
#endif
    
}

void
BamFormat::writeBamHeader()
{
#ifdef HAVE_HTSLIB
  int hdrWriteRes = sam_hdr_write(hFile, head);

#else
  
#endif
}
