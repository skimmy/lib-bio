#ifndef KMER_SCORE_H
#define KMER_SCORE_H

#include <vector>

#include "types.hpp"

/**
   \file kmerscore.hpp
   \brief This file contains function to compute the kmer score as defined in the
   thesis.

   The file also contains auxiliary types and constant definitions used throughout
   the `kmerscore` and (more generally) `algorithms` modules of the library.
 */


/**
   \fn isKmerUnique(const ReadkmerMap& map);
 */
bool isKmerUnique(const KmersMap& map);


#endif
