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
   \fn isKmerUniquelyMapped(const KmersMap& map);
   \brief Checks wether the passed mapping defines a unique position.

   This function can be used to discard k mers that have non unique mapping, the
   test checks whether all the aligned positions (<i>i.e.</i> those with a value
   different from `-1`) are consistent with each other and define a unique map
   position for the corresponding read.

   More precisely let `i` be the firs position with value different from `-1`,
   then all positions `j` of the input map must satisfy:
   `(map[j] == map[i] + (j - i))`

   Note that, since we ignore all the positions with mapping `-1`, a map with
   all such values (<i>i.e.</i> `map[i] = -1` for all `i`) pass the test and
   this function would return `true`. To keep consistency with this (counter
   intuitive behavior) the functions returns `true` whether an empty `map` is
   passed as its input.
 */
bool isKmerUniquelyMapped(const KmersMap& map);

std::vector< uint64_t > kmerScoreVector(const KmersMap& map, size_t k);

KmerScoreType scoreForVector(const std::vector< uint64_t >& v, size_t k);

size_t kmerErrorCount(const KmersMap& map, size_t k);

#endif
