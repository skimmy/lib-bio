// lcs.hpp

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

#ifndef SIM_LCS_H
#define SIM_LCS_H

#include <include/common.hpp>

#include <vector>
#include <iterator>

namespace lbio { namespace sim { namespace lcs {

template <typename _IterT>
std::vector<std::pair<lbio_size_t, lbio_size_t> > 
lcs_greedy(_IterT b1, _IterT e1, lbio_size_t n,
	   _IterT b2, _IterT e2, lbio_size_t m) {
  std::vector<std::pair<lbio_size_t, lbio_size_t>> lcs_list {};
  // TODO: NEEDS DEBUG
  lbio_size_t i {0};
  lbio_size_t j {0};
  //  while(b1 != e1 and b2 != e2) {
  while(i < n and j < m) {    
    if (*b1 == *b2) {
      lcs_list.push_back(std::make_pair(i, j));    
      ++i; ++b1;
      ++j; ++b2;
    }
    else {
      lbio_size_t d = 1;
      bool end_loop = false;
      while ( (!end_loop) and (d < (n-i+1) + (m-j+1)) ) {
	for (lbio_size_t h = 0; h < d; ++h) {
	  lbio_size_t di = d - h;
	  lbio_size_t dj = h;
	  if (*(b1+di) == *(b2+dj)) {
	    i += di;
	    j += dj;
	    lcs_list.push_back(std::make_pair(i,j));
	    b1 += di;
	    b2 += dj;
	    end_loop = true;
	  }
	}
	++d;
      }

    }
    
  }
  return lcs_list;
}

template <typename _IterT>
std::vector<std::pair<lbio_size_t, lbio_size_t> > 
lcs_greedy(_IterT b1, _IterT e1, _IterT b2, _IterT e2) {
  std::vector<std::pair<lbio_size_t, lbio_size_t>> lcs_list {};
  lbio_size_t n = std::distance(b1, e1);
  lbio_size_t m = std::distance(b2, e2);
  return lcs_greedy(b1, e1, n, b2, e2, m);
}

} } } 

#endif
