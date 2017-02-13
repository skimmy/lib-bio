// matrix.hpp
// Matrix classes and utilities

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

/**
   \file strctures/matrix.hpp
   Contains classes, structures and utilities to manipulate matrices.
 */

#ifndef LBIO_MATRIX_HPP
#define LBIO_MATRIX_HPP

#include <lbio.h>

DEFAULT_NAMESPACE_BEGIN

/**
   \brief A \c struct representing two dimensional matrix of generic type.

   \tparam _ContentT  The type of the elements stored in the matrix
   \tparam _Alloc     An allocator type

   This implementation uses one single array to store the two
   dimensional matrix. This means that accesses all happen in constant
   time, but moving rows and columns always requires to move the
   actual elements.
 */
template <typename _ContentT, typename _Alloc = std::allocator<_ContentT> >
struct _2D_matrix
{

  typedef _ContentT              content_type;
  typedef _ContentT*             content_pointer;
  typedef size_t                 size_type;
  typedef size_t                 index_type;
  
  // member variables
  content_pointer _mat;
  size_type _rows;
  size_type _cols;

  
  _2D_matrix(size_type _r, size_type _c)
    : _rows {_r}, _cols {_c}
  {
    _mat = new content_type[_rows*_cols];
  }

  ~_2D_matrix() {
    delete[] _mat;
  }

  content_type operator()(index_type _i, index_type _j) const {
    return _mat[_i*_cols + _j];
  }

  content_type& operator()(index_type _i, index_type _j) {
    return _mat[_i*_cols + _j];
  }
};

DEFAULT_NAMESPACE_END

#endif 

