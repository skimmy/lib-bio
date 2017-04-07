// generator.hpp

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

#ifndef LBIO_GENERATOR_HPP
#define LBIO_GENERATOR_HPP

#include <lbio.h>

#include <random>

DEFAULT_NAMESPACE_BEGIN

/**
   \brief This class implements an <em>IID sampler</em> based on the
   discrete probability distribution that is used as template \c _Dist

   \tparams The DiscreteProbability (or equivalent) class describing
   the sampled distribution
 */
template <typename _Dist>
class IIDSampler
{
public:
  using distribution = _Dist;
  using distribution_ref = distribution&;
  using sample_type = typename distribution::event_type;

  IIDSampler(distribution_ref d)
    : _std_distr {d.to_std_discrete_distribution()}  { }

  sample_type sample() {
    sample_type s = sample_type{0};
    return s;
  }
  
private:
  std::discrete_distribution<> _std_distr;
  
};

DEFAULT_NAMESPACE_END

#endif
