#include <iostream>
#include <iterator>

#include "../io.h"

void logInfo(const std::string& msg) {
  std::cerr << "[Info] " << msg  << std::endl;
}

void testFastqRead(const std::string& fastqPath) {
  logInfo("Fastq test on file " + fastqPath);

  std::cerr << "Istream iterator test...\n";
  std::ifstream ifs("/tmp/in.fastq", std::ifstream::in);

  for (std::istream_iterator<FastqRead> it(ifs);
       it != std::istream_iterator<FastqRead>(); ++it) {
    std::cout << *it;
  }

  ifs.close();
}

int main(int argc, char** argv) {

  testFastqRead(argc > 1 ? std::string {argv[1]} :
		std::string {"/tmp/in.fastq"} );
    
  return 0;
}
