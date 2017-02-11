#ifndef LBIO_MATRIX_HPP
#define LBIO_MATRIX_HPP

#include <lbio.h>

DEFAULT_NAMESPACE_BEGIN

template <typename _ContentT>
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

