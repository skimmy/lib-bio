// generator_task.hpp

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

#ifndef BIOTK_GENERATE_TASK_HPP
#define BIOTK_GENERATE_TASK_HPP

#include <lbio.h>
#include <prob/generator.hpp>

template <typename _Sampler>
std::string
generate_iid_bases(lbio_size_t len, _Sampler& sampler) {
  std::string _bases(len, 'N');
  sampler.sample(_bases.begin(), _bases.end());
  return _bases;
}

/**
   \brief Creates a new sequence of type _ContT from the input one by
   performing edit transformation according to array p which should
   indicate probability of substitution, deletion and insertion (match
   is obtained by complement) and, for substitution and insertions,
   output symbols are chosen i.i.d. among possible values.
 */
template <typename _ContT, typename _ArrayT>
_ContT
generate_iid_edit_transformation(_ContT _in, _ArrayT p) {
  
  // operation distribution
  std::discrete_distribution<int> _op_dist({p[0], p[1], p[2], 1-(p[0]+p[1]+p[2])});
  // iid for sub
  std::discrete_distribution<int> _sub_dist({1,1,1});
  // iid for ins
  std::discrete_distribution<int> _ins_dist({1,1,1,1});
  auto& gen = lbio::global_rand_generator<std::mt19937>();

  const std::vector<char> _bases {'A', 'C', 'G', 'T'};
  std::map<char,int> _rev_bases
  { {'A', 0}, {'C', 1}, {'G', 2}, {'T',3} };
  
  _ContT out {};
  auto _back_ins = std::back_inserter(out);
  auto _in_it = _in.cbegin();
  while(_in_it != _in.cend()) {
    int op = _op_dist(gen);
    std::cout << op;
    switch (op)
      {
      case 0: // Sub
	*_back_ins = _bases[ (_rev_bases[*_in_it] + _sub_dist(gen)) % 4];
	++_in_it;
	break;
	
      case 1: // Del
	++_in_it;
	break;
      case 2: // Ins
	*_back_ins = _bases[_ins_dist(gen)];
	break;
      case 3: // Match
	*_back_ins = *_in_it;
	++_in_it;
	break;
      default: // should never reach
	*_back_ins = *_in_it;
	++_in_it;
	break;
      }

  }
  return out;
}

#endif
