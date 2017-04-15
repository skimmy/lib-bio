// hamming_distance.hpp

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

#ifndef LBIO_HAMMING_DISTANCE_HPP
#define LBIO_HAMMING_DISTANCE_HPP

#include <lbio.h>

DEFAULT_NAMESPACE_BEGIN

template <typename _IterT>
lbio_size_t
hamming_distance(_IterT b1, _IterT e1, _IterT b2) {
  lbio_size_t d = 0;
  for (; b1 != e1; ++b1, ++b2) {
    d += static_cast<lbio_size_t>(*b1 != *b2);
  }
  return d;
}

DEFAULT_NAMESPACE_END

#endif

