#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <unordered_map>
#include <list>

#include "../sequence.h"

/**
   \file spectrum.hpp
   \brief Contains several helper functions to work with <i>spectrum</i> of a
   sequence.
 */


/**
   Defines a new type `KmerMap` to indicate as an `unordered_map` with key 
   `uint64_t`(for numeric k-mer) and value a list of positions (<i>i.e.</i>
   `size_t`)
 */
typedef std::unordered_map< uint64_t, std::list< size_t > > KmerMap;

/**
   \fn spectrumAsIntMap(const Sequence& ref, size_t k);
   \brief Computes the \f$ k \f$ spectrum of the given Sequence.

   The function uses the numeric version of a k-mer where the sequence of k
   nucleotides is represented as a sequence of \f$ 2k \f$ bits where each pair
   of bits represents one single nucleotide.

   Since <i>k-mers</i> are represented using a `uint64_t` value and two bits are
   used to represent a single nucleotide, this function does not allow k-mers
   with \f$ k > 32 \f$.

   The algorithms linearly scans the k-mers of the sequence a keep tracks of the
   amount of times each k-mer was observed. Since it uses an `urdoreed_map` to
   keep such index, the overall asymptotic complexity of this function is   
   \f[
   T = O(N T_{map})
   \f]
   where \f$ N \f$ is the length of the input sequence and \f$ T_{map} \f$ is
   the time needed to access (<i>i.e.</i> retrieve, modify and store) an element
   from the `unordered_map`.
   
   \param ref The reference as a Sequence type on which compute the spectrum
   \param k The size (<i>i.e.</i> number of symbols) of the k-mers
   
   \return a `unordered_map` from the STL containing [kmer, count] pairs where,
   for each kmer contained in ref, there is the total number of occurances of
   such kmer in the whole sequence.
 */
std::unordered_map< uint64_t, uint64_t > spectrumAsIntMap(const Sequence& ref, size_t k);

/**
   \fn kmersMapping(const Sequence& ref, size_t k);
   \brief Computes the mapping for all k-mers in a Sequence

   The function uses the numeric version of a k-mer where the sequence of k
   nucleotides is represented as a sequence of \f$ 2k \f$ bits where each pair
   of bits represents one single nucleotide.

   Since <i>k-mers</i> are represented using a `uint64_t` value and two bits are
   used to represent a single nucleotide, this function does not allow k-mers
   with \f$ k > 32 \f$.

   The algorithms linearly scans the k-mers of the sequence a keep tracks of the
   positions where each k-mer was observed. Since it uses an `urdoreed_map` to
   keep such index, the overall asymptotic complexity of this function is   
   \f[
   T = O(N T_{map} T_{vec})
   \f]
   where \f$ N \f$ is the length of the input sequence, \f$ T_{map} \f$ is
   the time needed to access (<i>i.e.</i> retrieve, modify and store) an element
   from the `unordered_map` and \f$ T_{vec} \f$ is the time needed to add an
   element to a `vector`.

   \param ref The reference as a Sequence type on which perform the mapping
   \param k The size (<i>i.e.</i> number of symbols) of the k-mers

   \return a `unordered_map` containing, for each of the k-mer in the reference sequence,
   a `vector` of positions where that k-mer was observed

 */
KmerMap kmersMapping(const Sequence& ref, size_t k);



#endif
