#include <include/common.hpp>

#include <include/options.hpp>

#include <include/align.hpp>
#include <include/io.hpp>
#include <include/log.hpp>

#include <list>
#include <iostream>

void
evaluateAlignmentScore(const Options &opts) {
    if (opts.inputReference.empty() || opts.inputSAM.empty()) {
        fatal_error("(AlignScore): invalid reference and/or SAM");
    }

    // load reference
    std::string ref = loadFromFile(opts.inputReference);
    std::cout << ref.size() << std::endl;

    // load SAM
    std::list<AlignPair> aligns;
    loadAlignFromSAM(opts.inputSAM, aligns);
    std::cout << aligns.size() << std::endl;
}
