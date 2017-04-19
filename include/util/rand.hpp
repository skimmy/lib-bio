// rand.hpp

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

#ifndef LBIO_RAND_HPP
#define LBIO_RAND_HPP

#include <lbio.h>

#include <random>
#include <chrono>


DEFAULT_NAMESPACE_BEGIN

template <typename _GenT>
_GenT& global_rand_generator() {
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  static _GenT g1 (seed);
  return g1;

}

template <typename _GenT>
_GenT& global_rand_generator_seed(unsigned int seed) {
  static _GenT g1 (seed);
  return g1;

}


DEFAULT_NAMESPACE_END

#endif
