#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <unordered_map>

#include "../sequence.h"

/**
   \file spectrum.hpp
   \brief Contains several helper functions to work with <i>spectrum</i> of a
   sequence.
 */


/**
   \fn spectrumAsIntMap(const Sequence& ref);
   \brief Computes the \f$ k \f$ spectrum of the given Sequence.

   The function uses the numeric version of a k-mer where the sequence of k
   nucleotides is represented as a sequence of \f$ 2k \f$ bits where each pair
   of bits represents one single nucleotide.

   Since <i>k-mers</i> are represented using a `uint64_t` value and two bits are
   used to represent a single nucleotide, this function does not allow k-mers
   with \f$ k > 32 \f$.
 */
std::unordered_map< uint64_t, uint64_t > spectrumAsIntMap(const Sequence& ref);

#endif
