#include "FastFormat.hpp"
#include "FastqFormat.hpp"

#include <../core/Reference.hpp>

#include <string>
#include <iostream>

int main(int argc, char** argv) {
  std::string dnaPath("/storage/bio/ecoli/ecoli.fasta");
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
  std::cout << reads.getNextRead() << std::endl;
  return 0;
}
