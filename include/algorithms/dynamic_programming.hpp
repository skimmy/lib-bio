#ifndef DYNAMIC_PROGRAMMING_HPP_
#define DYNAMIC_PROGRAMMING_HPP_

#define DEFAULT_NAMESPACE_BEGIN namespace lbio {
#define DEFAULT_NAMESPACE_END }  

DEFAULT_NAMESPACE_BEGIN

template <typename _ContentT>
struct dp_matrix
{
  // Convenienve type re-definition (inspired by the std container impl)
  typedef _ContentT      content_type;
  typedef content_type*  content_pointer;
  typedef size_t         size_type;
  typedef size_t         index_type;

  // members 
  content_pointer _mat;
  size_type _rows;
  size_type _cols;

  // constructor accepting dimensions of the matrix
  dp_matrix(size_type _r, size_type _c)
    : _rows {_r}, _cols {_c}
  {
    // TODO: Eventually convert to allocator template
    _mat = new content_type[_rows*_cols];
  }

  ~dp_matrix() {
    delete[] _mat;
  }

  // elements are accessed via M(i,j) rather than M[i][j]
  content_type& operator()(index_type i, index_type j) {
    return _mat[i*_cols + j];
  }
  content_type operator()(index_type i, index_type j) const {
    return _mat[i*_cols + j];
  }
      

};

DEFAULT_NAMESPACE_END

#endif
