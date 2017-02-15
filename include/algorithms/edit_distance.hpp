// edit_distance.hpp

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

#ifndef LBIO_EDIT_DISTANCE_HPP
#define LBIO_EDIT_DISTANCE_HPP

#include <lbio.h>
#include <algorithms/dynamic_programming.hpp>

#include <iterator>

DEFAULT_NAMESPACE_BEGIN

/**
   \brief Implements the Wagner Fischer algorithm for the edit distance.

   \tparam _SeqT  The sequence type on which find the edit distnace
   \tparam _CostT The cost type 
 */
template<typename _SeqT, typename _CostT = lbio_size_t>
class edit_distance_wf
{
public:
  typedef _SeqT                  sequence_type;
  typedef _CostT                 cost_type;
  typedef dp_matrix<cost_type>   matrix;
  typedef lbio_size_t            size_type;
  typedef lbio_size_t            index_type;

  edit_distance_wf(size_type _n, size_type _m)
    : _mat {_n+1, _m+1} {
    init();
  }
  
  cost_type
  compute(const sequence_type& _seq1, size_type _n1,
	  const sequence_type& _seq2, size_type _n2) {
    auto _iter1 = _seq1.begin();
    for (index_type _i=1; _i<=_n1; ++_i, ++_iter1) {
      auto _iter2 = _seq2.begin();
      for (index_type _j=1; _j<=_n2; ++_j, ++_iter2)  {
	cost_type delta = (*_iter1 == *_iter2) ? 0 : 1;
	_mat(_i, _j) = std::min(_mat(_i-1,_j-1) + delta,
				std::min(_mat(_i-1,_j) + 1, _mat(_i,_j-1) + 1));
      }
    }
    
    return _mat(_n1, _n2);
  }

  cost_type
  compute(const sequence_type& _seq1, const sequence_type& _seq2) {
    size_type n = std::distance(_seq1.begin(), _seq1.end());
    size_type m = std::distance(_seq2.begin(), _seq2.end());
    return compute(_seq1, n, _seq2, m);
  }

private:
  matrix _mat;

  void
  init() {
    for (index_type _i=0; _i<_mat._rows; ++_i) {
      _mat(_i,0) = _i;
    }
    for (index_type _j=0; _j<_mat._cols; ++_j) {
      _mat(0,_j) = _j;
    }
  }
};

/**
   \brief This struct represents an edit transformation as a pair x ->
   y of the generic type.
   
   \tparam _TransfT This is the type used to represent transformations

   \note The type used to instantiate the template must support
   assignment and equality comparison. Moreover it must be
   initializable from the integral value 0 (used in the default
   constructor).
   
 */
template<typename _TransfT>
struct edit_transformation
{
  typedef _TransfT    transformation;

  const transformation _epsilon;

  transformation _from;
  transformation _to;

  edit_transformation()
    : _epsilon {0}, _from{_epsilon}, _to{_epsilon} { } 

  explicit edit_transformation(transformation _f, transformation _t,
			       transformation _e = 0)
    : _epsilon {_e}, _from {_f}, _to {_t} { }

  void
  assign(transformation _f, transformation _t) {
    _from = _f;
    _to = _t;
  }

  void
  checked_assign(transformation _f, transformation _t) {
    if ( (_f == _epsilon) && (_t == _epsilon) ) {
      throw domain_error("Invalid edit transformation");
    }
    assign(_f,_t);
  }

  bool
  is_match() const {
    return ( (_from == _to) && (_to != _epsilon) );
  }

  bool
  is_substitute() const {
    return ( (_from != _to) && (_from != _epsilon) && (_to != _epsilon) );
  }

  bool
  is_delete() const {
    return ( (_from != _epsilon) && (_to == _epsilon) );
  }

  bool
  is_insert() const {
    return ( (_from != _epsilon) && (_to == _epsilon)  );
  }

  bool
  is_valid() const {
    return ( (_from != _epsilon) && (_to != _epsilon) );
  }
 
};

DEFAULT_NAMESPACE_END

#endif
