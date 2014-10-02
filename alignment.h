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
#include "alignment/aligndef.hpp"

std::vector<Position<int>> alignReads(const std::vector<Read>& reads, const std::string& ref);


#endif
