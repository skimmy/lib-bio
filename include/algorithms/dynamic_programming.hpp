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
