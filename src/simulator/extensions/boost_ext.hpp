#ifndef SIM_EXT_BOOST_H
#define SIM_EXT_BOOST_H

// If autotools tells us that boost is available we define extnsions
#if defined(HAVE_BOOST)

// This can be used to check if the boost extensions are avialble
#define SIM_BOOST_EXTENSIONS

// Each extension goes into a proper include for better modularity

#include "boost_ext/boost_edit.hpp"

#endif

#endif
