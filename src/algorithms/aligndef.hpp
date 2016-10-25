#ifndef ALIGN_DEF_H
#define ALIGN_DEF_H

#include "../core.h"
#include "Position.hpp"

#include <vector>
#include <string>


enum BacktrackOperation {Unset, Match, Substitution, Insertion, Deletion};

namespace bio {

std::vector< Position<int> > alignReads(const std::vector<Read>& reads, const std::string& ref);

}

#endif
