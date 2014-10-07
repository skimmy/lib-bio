#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <vector>
#include <string>

#include "sequence.h"

/**
 * All alignment inclusions here
 */
#include "alignment/SmithWatermanDP.hpp"
#include "alignment/Position.hpp"
#include "alignment/ScoredPosition.hpp"
#include "alignment/aligndef.hpp"

// Alignment functions
std::vector<Position<int>> alignReads(const std::vector<Read>& reads, const std::string& ref);

// K-Specrtum functions
std::vector<uint64_t>* getKmersFrequency(const Sequence& seq, size_t k);


#endif
