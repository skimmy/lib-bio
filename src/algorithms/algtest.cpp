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
  return 0;
}
