#ifndef TASKS_H
#define TASKS_H

#include <string>
#include <vector>

#include "alignment.h"

/** \file tasks.hpp
 * \brief Contains definition of functions used to perform different tasks based
 * on the input received from the command line.
 */

/**  \fn alignFastqReadsSimpleSW(const string& readsPath, const string& referencePath, uint nThreads = 1, size_t nReads = -1);
  \brief Aligns reads contained in a fastq file against a regerence contained in
  a fast file.

  The function opens the files at the paths indicated by input strings. 
  If any input operation fails, the function throws an ??? exception.
 
  This is a threaded function, caller can specify the number of threads that the
  algorithm is supposed to use, if nothing is specified the algorithm is run in
  aingle threaded mode.
 
  Please note that <i>this is a blocking function</i> that will not return until
  all threads have completed their execution.
 
  We call <i>reads block</i> all the reads assigned to the same thread, to decide
  decide the size of the blocks the function uses am heuristic that assumes each
  read contributes for a fixed amount of bytes B to the overall size of the file
  (which is reasonable for "standard" fastq files). In other words a fraction of
  \f$ 1/T \f$ of the total file is read before creating a new block \f$ T \f$ is
  the number of threads passed as parameter
 
  \param readsPath The path of the fastq file containing the reads
  \param referencePath The path to the fast file containing the reference sequence
  \param nThreads The number of threads used to align (can be omitted, default is 1)
  \param nReads The number of reads, can be omitted in which case this number is
  estimated when creating reads blocks (see above).
*/
std::vector<Position<int>> alignFastqReadsSimpleSW(const string& readsPath, const string& referencePath, uint64_t nThreads = 1, size_t nReads = -1);

#endif
