// io_helper.hpp

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

#ifndef LBIO_IO_HELPER_HPP
#define LBIO_IO_HELPER_HPP

#include <lbio.h>

#include <algorithm>
#include <iostream>
#include <map>


DEFAULT_NAMESPACE_BEGIN

template <typename _T>
_T from_string(std::string str) {
  std::stringstream _stream{ str };
  _T t;
  _stream >> t;
  return t;
}

template<>
inline std::string
from_string<std::string>(std::string str) {
  // hopefully compiler will optimize
  return  str;
}

template <typename _KeyT, typename _ValT>
std::map<_KeyT,_ValT>
stream_to_map(std::istream& _is, char delimiter=',', char comment='#') {
  std::map<_KeyT, _ValT> _map;
  while(!_is.eof()) {
    std::string line {};
    std::getline(_is, line);
    // empty lines and comment lines are allowed
    if (line.empty() || line[0]==comment) {
      continue;
    }
    // Tokenization: assumes 'Key [sep] Val' other occurrences of [sep]
    // are considered part of Val
    auto _sep_iter = std::find(line.begin(), line.end(), delimiter);
    _KeyT _key = from_string<_KeyT>(std::string {line.begin(), _sep_iter});
    _ValT _val = from_string<_ValT>(std::string { ++_sep_iter, line.end()});
    _map[_key] = _val;
  }
  return _map;
}

using StringStringMap = typename std::map<std::string, std::string>;

DEFAULT_NAMESPACE_END

#endif
