// Public library include str_util.hpp

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

#ifndef LBIO_STR_UTIL_HPP
#define LBIO_STR_UTIL_HPP

#include <lbio.h>

#include <vector>

DEFAULT_NAMESPACE_BEGIN

template <typename _It>
std::vector<std::string>
split(_It beg, _It end, char delim = ' ') {
  std::vector<std::string> v;
  _It s = beg;
  while(beg != end) {
    if ((*beg) == delim) {
      v.push_back(std::string(s, beg));
      s = beg;
      ++s;
    }
    ++beg;
  }
  if (s != beg) {
    v.push_back(std::string(s, beg));
  }
  return v;
}

DEFAULT_NAMESPACE_END

#endif
