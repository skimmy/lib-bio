#include "../algorithms.h"

#include <iostream>
#include <vector>

int main(int argc, char** argv) {
  Reference r1("ACGTT");
  Reference r2("AGGTT");
  std::cout << bio::hammingDistance(r1,r2) << std::endl;
  Read read1;
  read1.setBases("GTT");
  std::vector<Read> v(1);
  v[0] = read1;
  std::vector< Position<int> > p = bio::alignReads(v, "ACGTT");
  for (Position< int > pi : p) {
    std::cout << pi.getPosition();
  }
  std::cout << std::endl;
  uint64_t D2 = dstats::D2(r1, r2, 4);
  std::cout << D2 << std::endl;

  std::cout << "\n\n" << DNAAlphabet2Bits::charToInt('A') << "\n";
  return 0;
}
