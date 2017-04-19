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

#include <util/rand.hpp>

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
  using event_type = typename distribution::event_type;
  using probability_type = typename distribution::probability_type;

  /**
     \brief Constructs an IIDSampler from a \c distribution
   */
  IIDSampler(const distribution& _d)
    : _events {} {
    std::vector<probability_type> weights {};
    auto _d_map = _d.get_container();
    for (auto event : _d_map) {
      _events.push_back(event.first);
      weights.push_back(event.second);
    }
    _dist = std::discrete_distribution<>(weights.begin(), weights.end());
  }

  /**
     \brief Returns a \c std::vector of elements of type \c event_type
     containing \c k samples from the underlying distribution
   */
  std::vector<event_type>
  sample(lbio_size_t k = 1) {
    auto& _gen = global_rand_generator<std::mt19937>();
    std::vector<event_type> samples;
    while (k > 0) {
      samples.push_back(_events[_dist(_gen)]);
      --k;
    }
    return samples;
  }

  /**
     \brief Fills the range <em>[_beg, _end]</em> with random samples
     from the underlying distribution
   */
  template <typename _IterT>
  lbio_size_t
  sample(_IterT _beg, _IterT _end) {
    auto& _gen = global_rand_generator<std::mt19937>();
    lbio_size_t count {0};
    for (; _beg != _end; ++_beg, ++count) {
      *_beg = _events[_dist(_gen)];
    }
    return count;
  }

private:
  std::vector<event_type> _events;
  std::discrete_distribution<> _dist;  
};


DEFAULT_NAMESPACE_END

#endif
