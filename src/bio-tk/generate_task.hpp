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

#include <lbio.h>

template <typename _Sampler>
std::string
generate_iid_bases(lbio_size_t len, _Sampler sampler) {
  std::string _bases(len, 'N');
  sampler.sample(_bases.begin(), _bases.end());
  return _bases;
}
