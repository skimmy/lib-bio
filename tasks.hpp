#ifndef TASKS_H
#define TASKS_H

#include <string>
#include <vector>

#include "alignment.h"

/**
 * \brief Aligns reads contained in a fastq file against a regerence contained in
 * a fast file.
 *
 * \param readsPath The path of the fastq file containing the reads
 * \param referencePath The path to the fast file containing the reference sequence
 */
std::vector<Position<int>> alignFastqReadsSimpleSW(const string& readsPath, const string& referencePath);

#endif
