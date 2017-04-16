// d2_distance.hpp

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

#ifndef LBIO_D2_DISTANCE_HPP
#define LBIO_D2_DISTANCE_HPP

#include <lbio.h>

#include <unordered_map>
#include <limits>

DEFAULT_NAMESPACE_BEGIN

/**
   \brief Counts the number of k-mers in the range [_b,_e)
 */
template <typename SeqT, typename IterT, typename CountT = lbio_size_t>
std::unordered_map<SeqT, CountT>
kmer_count(lbio_size_t k, IterT _b, IterT _e) {
  std::unordered_map<SeqT, CountT> kmer_map;
  // first kmer
  IterT _l = _b;
  for (lbio_size_t kk = 0; kk < k; ++kk) {
    if (_b == _e) {
      return kmer_map; // empty map
    }
    ++_b;
  }
  kmer_map[SeqT(_l, _b)]++;
  while(_b != _e) {
    _b++;
    _l++;
    kmer_map[SeqT(_l, _b)]++;
  }
  return kmer_map;
}

/**
   \brief Computes \e D2 statistics on the ranges [b1,e1) and [b2,e2)
   using k as size of \e k-mer
 */
template <typename SeqT, typename IterT>
double
D2(lbio_size_t k, IterT b1, IterT e1, IterT b2, IterT e2) {
  double dist = 0;  
  // count kmers
  std::unordered_map<SeqT, lbio_size_t> kmer_count_1 =
    kmer_count<SeqT, IterT>(k, b1, e1);
  std::unordered_map<SeqT, lbio_size_t> kmer_count_2 =
    kmer_count<SeqT, IterT>(k, b2, e2);
  // 'dot product'
  for (auto x : kmer_count_1) {
    if (kmer_count_2.count(x.first)) {
      dist += x.second * kmer_count_2[x.first];
    }
  }
  return dist;
}

/**
   \brief Computes \e D2* statistics on the ranges [b1,e1) and [b2,e2)
   using k as size of \e k-mer.  To get null hypothesis probability
   uses the functional p.
 */
template <typename SeqT, typename IterT, typename ProbW>
double
D2_star(lbio_size_t k, IterT b1, IterT e1, IterT b2, IterT e2, ProbW p) {
  double dist = 0;
  // counts kmer
  std::unordered_map<SeqT,double> kmer_map_1 =
    kmer_count<SeqT, IterT, double>(k, b1, e1);
  std::unordered_map<SeqT,double> kmer_map_2 =
    kmer_count<SeqT, IterT, double>(k, b2, e2);
  // normalization of kmer count
  // Notice D2* assumes input sequences have same len n
  double barn_1 = std::distance(b1, e1) - k + 1;
  for (auto& x : kmer_map_1) {
    x.second = x.second - barn_1 * p(x.first);
  }
  for (auto& x : kmer_map_2) {
    x.second = x.second - barn_1 * p(x.first);
  }
  // dot product
  for (auto x : kmer_map_1) {
    if (kmer_map_2.count(x.first)) {
      dist += x.second * kmer_map_2[x.first]
	/ ( barn_1 * p(x.first) );

    }
  }
  return dist;
}

DEFAULT_NAMESPACE_END

#endif
