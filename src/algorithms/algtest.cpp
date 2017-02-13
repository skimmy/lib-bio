// algtest.cpp

// Copyright 2017 Michele Schimd

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//     http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#include "../algorithms.h"

#include <lbio.h>
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

  lbio_size_t n = 8;
  lbio_size_t m = 8;
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
