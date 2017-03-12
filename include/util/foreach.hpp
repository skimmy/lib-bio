// foreach.hpp

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

#ifndef LBIO_FOREACH_HPP
#define LBIO_FOREACH_HPP

#include <algorithm>

DEFAULT_NAMESPACE_BEGIN

template <typename _Iter1, typename _Iter2,
	  typename _Op, typename _Pred>
void inline
for_each_such_that(_Iter1 _beg,  _Iter2 _end,  _Pred _pred, _Op _oper)
{
  std::for_each(_beg, _end, [=](decltype(*_beg) val)
    { if (_pred(val)) {
	_oper(val);
      }
    });
}


DEFAULT_NAMESPACE_END

#endif
