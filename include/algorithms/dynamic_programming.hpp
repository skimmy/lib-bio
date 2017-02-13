// dynamic_programming.hpp

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

#ifndef DYNAMIC_PROGRAMMING_HPP_
#define DYNAMIC_PROGRAMMING_HPP_

#include <lbio.h>
#include <structures/matrix.hpp>

DEFAULT_NAMESPACE_BEGIN

template <typename _ContentT>
struct dp_matrix : public _2D_matrix<_ContentT>
{
  // Convenience type re-definition (inspired by the std container impl)
  typedef _2D_matrix<_ContentT>  _Base;
  typedef _ContentT              content_type;
  typedef _ContentT*             content_pointer;
  typedef size_t                 size_type;
  typedef size_t                 index_type;

  
  dp_matrix(size_type _r, size_type _c)
    : _Base(_r,_c)
  {
  }
};

DEFAULT_NAMESPACE_END

#endif
