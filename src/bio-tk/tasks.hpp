#ifndef TASKS_H
#define TASKS_H

#include <cstdint>

#include <string>
#include <vector>

#include "../algorithms.h"

/** \file tasks.hpp
 * \brief Contains definition of functions used to perform different tasks based
 * on the input received from the command line.
 */

/**  \fn alignFastqReadsSimpleSW(const string& readsPath, const string& referencePath, std::ostream& output, uint64_t nThreads = 1, size_t nReads = -1);
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
std::vector<ScoredPosition<int,int> > alignFastqReadsSimpleSW(const string& readsPath, const string& referencePath, std::ostream& output, uint64_t nThreads = 1, size_t nReads = -1);

/**
   \fn taskComputeKSpectrum(size_t k, const string& referenceFile);
   \brief Computes the \f$ k \f$ spectrum of the input reference file.
 */
void taskComputeKSpectrum(size_t k, const string& referenceFile);

/**
   \fn taskMapReadsKmers(const string& reference, const string& reads, size_t k, const string& out = "", size_t nThreads = 1);
   \brief For each read in the input set, maps its kmers against the given reference.

   This function takes as input file name for reference sequence and reads set, the size
   of kmers \f$ k \$f and a (possibly empty) `string` representing the output path.
   First it creates a _map_g
 */
void taskMapReadsKmers(const string& reference, const string& reads, size_t k, const string& out = "");

void taskKmerScoreReads(const string& reference, const string& reads, size_t k, const string& out = "", size_t nThreads = 1);

void task_read_statistics(const std::string& reads, const string& w_dir, const string& prefix);

void
task_generate(std::map<std::string,std::string> gen_params);


#endif
