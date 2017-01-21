#ifndef SIM_BOOST_EXT_EDIT_H
#define SIM_BOOST_EXT_EDIT_H

// This file contains the extensions for edit distance module that
// uses boost features: uBlas, ...

// The following guard is set in common.hpp therefore when defined we
// can assume that common.hpp has been included and that basic
// definition have already happened (e.g., lbio_size_t)
#if defined(LBIO_EDIT_BOOST_EXTENSIONS)

#include <boost/numeric/ublas/matrix.hpp>

#define ublas_ boost::numeric::ublas

template<typename CostType>
struct DynamicProgrammingBoost {
  lbio_size_t n;
  lbio_size_t m;
  ublas_::matrix<CostType> dp_matrix;


  DynamicProgrammingBoost(lbio_size_t n_, lbio_size_t m_)
    : n(n_), m(m_), dp_matrix(n_,m_) { }
};

#endif

#endif
