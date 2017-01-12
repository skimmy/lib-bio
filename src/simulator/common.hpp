#ifndef COMMON_H
#define COMMON_H

// This is config.h generated using autotools. It contains all the variables
// that should be set at any time and therefore it is at the top of this
// header rather than in the inclusion section.
#include <config.h>


///////////////////////////////////////////////////////////////////////////
//
//           DEFINITION OF SOME USEFUL MACROS AND CONSTANTS
//
///////////////////////////////////////////////////////////////////////////

// Uncomment if NDEBUG is needed and not generated during compiler
//#define NDEBUG 


// Define a NOP_MACRO to be used to expand no operation. For example this
// is used to mask output when NDEBUG is set.
// NB. Code is inspired by <assert.h> (GNU related guard has been removed)
#if defined __cplusplus 
# define LBIO_NOOP static_cast<void>
#else
# define LBIO_NOOP (void)
#endif



///////////////////////////////////////////////////////////////////////////
//
//                  HERE ARE ALL THE INCLUSIONS NEEDED
//
///////////////////////////////////////////////////////////////////////////

// It is bad design to use a single file to include every header in the
// project. Such solution is, however, convenient for small projects.


// includes for simulator
#include "generator.hpp"
#include "online.hpp"
#include "options.hpp"
#include "chain.hpp"
#include "prob.hpp"
#include "util.hpp"
#include "align.hpp"
#include "io.hpp"
#include "edit.hpp"
#include "log.hpp"

// this has likely been included in some of the above files. We just keep
// it because it is explictly required by definition of function(s) below
#include <string>

// testing functions
void testAll();
void fatal_error(const std::string &msg, int exit_code = 1);

#endif
