#include "../algorithms.h"

#include <algorithms/dynamic_programming.hpp>

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

  size_t n = 8;
  size_t m = 8;
  lbio::dp_matrix<int> matrix {n,m};
  for (int i = 0; i < matrix._rows; ++i) {
    for (int j = 0; j < matrix._cols; ++j) {
      matrix(i,j) = i - j ;
    }
  }
  for (int i = 0; i < matrix._rows; ++i) {
    for (int j = 0; j < matrix._cols; ++j) {
      std::cout << matrix(i,j) << "\t" ;
    }
    std::cout << "\n";
  }
  return 0;
}
