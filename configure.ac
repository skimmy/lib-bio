AC_INIT([libbio], [1.0], [])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CXX
AC_PROG_RANLIB
AM_PROG_AR

AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([memset])
AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([sqrt])
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT64_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_TYPE_SIZE_T


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
  Makefile
  src/core/Makefile
  src/io/Makefile
  src/algorithms/Makefile
  src/simulator/Makefile
  src/bio-tk/Makefile
])
AC_OUTPUT 
