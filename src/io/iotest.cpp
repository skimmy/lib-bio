#include "BamFormat.hpp"
#include <string>
#include <iostream>
#include <memory>

typedef std::unique_ptr<UIntList> pUIntList;

int main(int argc, char** argv) {
  /*  std::string dnaPath("/storage/bio/ecoli/ecoli.fasta");
  FastFormat ecoli(dnaPath);
  Reference r = ecoli.toReference();
  std::cout << std::endl;
  for (int i = 0; i < 102; ++i) {
    //    std::cout << r.getBaseAt(i);
  }
  std::cout << std::endl;
  std::string readsPath("/storage/bio/ecoli/ecoli_sample.fastq");
  FastqFormat reads;
  reads.openFile(readsPath);
  std::cout << reads.getNextRead() << std::endl;
  std::cout << reads.getNextRead() << std::endl;
  std::cout << reads.getNextRead() << std::endl;
  std::cout << reads.getNextRead() << std::endl;*/

  lbiobam::BamFormat bam;
  bam.open("/tmp/sample.bam");
  IdPos idPos = bam.getNext();
  while (idPos.first != "EOF") {
    std::cout << idPos.first << " " << idPos.second << "\n";
    idPos = bam.getNext();
  }
  /*  pUIntList pList = bam.getAlignmentPositions();
  std::cout << "Found " << pList->size() << " alignments\n";
  for (uint64_t align : *pList) {
    std::cout << align << "\n";
    }*/

  // BAM write tests
  lbiobam::BamFormat oBam;
  oBam.open("/tmp/out.sam", lbiobam::BamOpenWrite);
  // ... DO STUFF ...
  oBam.close();
  
  std::cout << std::endl;
  bam.close();
  
  return 0;
}
