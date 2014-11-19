#ifndef ALGS_TYPES_H
#define ALGS_TYPES_H

#include <cstdint>

#include <vector>
#include <list>
#include <unordered_map>

/**
   \file types.hpp
   \brief Defines several types for different part of the algorithm module.
 */


/**
   \brief Defines a type for mapping ordered k-mers into <b>one</b> single
   position.

   The mapping of k-mers in a sequence of length \f$ m \f$ is a vector of length
   \f$ \bar{m} = m - k + 1 \f$, element at index  \f$ i \in [0, \bar{m}[ \f$ contains
   either a position on a reference sequence (<i>i.e.</i> in \f$ [0, n[ \f$ ) or `-1`
   to indicate no map for the corresponding position.

   The fact that mapped positions are implicitly ordered on a `vector` structure
   can be used to map consecutive kmers of the input sequence. 
 */
typedef std::vector< size_t > KmersMap;

/**
   \brief The k-mer score type defined as a long (may change in the future)
 */
typedef long KmerScoreType;

/**
   \brief Defines a type for indexing numeric k-mers in a sequence.
   
   Numeric k-mers are here represented in their `uint64_t` form so to be suitable
   for key of maps, the value part is a `std::list` of `size_t` elements that are
   used to indicate the position(s) where the k-mer maps.

 */
typedef std::unordered_map< uint64_t, std::list< size_t > > NumericKmerIndex;

/**
   \brief Map taking `uint64_t` as both key and value.
 */
typedef std::unordered_map< uint64_t, uint64_t > IntMap;


#endif
