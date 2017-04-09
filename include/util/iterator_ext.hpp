// probabilitiesx.hpp

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

#ifndef LBIO_ITERATOR_EXT__HPP
#define LBIO_ITERATIR_EXT_HPP

#include <lbio.h>


DEFAULT_NAMESPACE_BEGIN

/**
   \brief Iterates thorugh elements of a \c map, dereferencing only
   keys

   \tparams _Map the underliying \c map
 */
template <typename _Map>
class map_key_iterator
{
public:
  // this iterators have deference key of the underlying map
  using value_type = typename _Map::key_type;
  using map_iterator = typename _Map::iterator;
  
  map_key_iterator(_Map& _orig_map) {
    _map_iter = _orig_map.begin();
  }

  value_type
  operator*() {
    return (*_map_iter).first;
  }
  
  map_key_iterator<_Map>&
  operator++() {
    ++_map_iter;
    return *this;
  }

  bool
  operator != (const map_key_iterator<_Map>& other) const {
    return (_map_iter != other._map_iter);
  }

  bool
  operator != (const map_iterator& _map_it) const {
    return _map_iter != _map_it;
  }
  
private:
  typename _Map::iterator _map_iter;
};

/**
   \brief Iterates thorugh elements of a \c map, dereferencing only
   values

   \tparams _Map the underliying \c map
 */
template <typename _Map>
class map_value_iterator
{
  public:
  // this iterators have deference key of the underlying map
  using value_type = typename _Map::mapped_type;
  using map_iterator = typename _Map::iterator;
  
  map_value_iterator(_Map& _orig_map) {
    _map_iter = _orig_map.begin();
  }

  value_type
  operator*() {
    return (*_map_iter).second;
  }
  
  map_value_iterator<_Map>&
  operator++() {
    ++_map_iter;
    return *this;
  }

  bool
  operator != (const map_value_iterator<_Map>& other) const {
    return (_map_iter != other._map_iter);
  }

  bool
  operator != (const map_iterator& _map_it) const {
    return _map_iter != _map_it;
  }
  
private:
  map_iterator  _map_iter;
};


DEFAULT_NAMESPACE_END

#endif
