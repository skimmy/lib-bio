#include <include/common.hpp>
#include <include/online.hpp>

#include <include/generator.hpp>


#include <cstring>

using namespace lbio::sim::generator;

GenomeSegment::GenomeSegment(size_t N, size_t m, size_t l) {
    this->length = l;
    this->total_length = N;
    this->read_length = m;
    this->current_start = 0;

    if (N < l) {
        this->length = N;
    }

    this->genome = new char[length];

}

GenomeSegment::~GenomeSegment() {
    if (this->genome != NULL) {
        delete[] this->genome;
    }
}

void generateFirstGenomeSegment(GenomeSegment &g) {
    generateIIDGenome(g.length, g.genome);
}

void generateNewGenomeSegment(GenomeSegment &g, size_t keep_pref) {
    // copy the last keep_pref character in the first one (to ensure pending reads are
    // corrected generated also after regeneration of the genome
    strncpy(g.genome, g.genome + (g.length - keep_pref), keep_pref);

    // generate new genome
    generateIIDGenome(g.length - keep_pref, g.genome + keep_pref);
}
