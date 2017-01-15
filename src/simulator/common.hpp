#ifndef COMMON_H
#define COMMON_H

// This is config.h generated using autotools. It contains all the
// variables that should be set at any time and therefore it is at the
// top of this header rather than in the inclusion section.
#include <config.h>

// this is needed for definition of functions below and also for the
// definition of some types (e.g., size_t) which are used for some
// type renaming/using/expansion (e.g., lbio_size_t)
#include <string>


//////////////////////////////////////////////////////////////////////
//
//         DEFINITION OF SOME USEFUL MACROS AND CONSTANTS
//
//////////////////////////////////////////////////////////////////////

// Uncomment if NDEBUG is needed and not generated during compiler
//#define NDEBUG 


// Define a NOP_MACRO to be used to expand no operation. For example
// this is used to mask output when NDEBUG is set.  NB. Code is
// inspired by <assert.h> (GNU related guard has been removed)
#if defined __cplusplus 
# define LBIO_NOOP static_cast<void>
#else
# define LBIO_NOOP (void)
#endif

// Define a custom type for size type, this is just a convenience and
// it will take long time (if ever) to port every part of the
// code. For the time being this type exapnds to size_t, but in the
// future changes there may be needed.
#if defined __cplusplus 
using lbio_size_t = size_t;
#else
#define lbio_size_t size_t
#endif

//////////////////////////////////////////////////////////////////////
//
//                HERE ARE ALL THE INCLUSIONS NEEDED
//
//////////////////////////////////////////////////////////////////////

// It is bad design to use a single file to include every header in 
// the project. Such solution is, however, convenient for small
// projects.


// includes for simulator
// #include "generator.hpp"
// #include "online.hpp"
// #include "options.hpp"
// #include "chain.hpp"
// #include "prob.hpp"
// #include "util.hpp"
// #include "align.hpp"
// #include "io.hpp"
// #include "edit.hpp"
// #include "log.hpp"

// The followinf arrays are used to translate int to bases and
// vice-versa. They are declare constexpr to allow compiler
// optimizations (when possible).
constexpr char bases[] = {'A', 'C', 'G', 'T'};
// all but ACGTacgt are mapped to 127
constexpr char revBases[] = {
  127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 
  127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 
  127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 
  127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 
  127, 0,   127, 1,   127, 127, 127, 2,   127, 127, 127, 127, 127, 127, 127, 127, 
  127, 127, 127, 127, 3,   127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 
  127, 0,   127, 1,   127, 127, 127, 2,   127, 127, 127, 127, 127, 127, 127, 127, 
  127, 127, 127, 127, 3, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
};



// testing functions
void testAll();
void fatal_error(const std::string &msg, int exit_code = 1);

#endif
