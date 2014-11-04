#ifndef OPTIONS_H
#define OPTIONS_H

#include <iostream>

using namespace std;

enum GenomeFormat { GENOME_CUSTOM, FAST };
enum ReadsFormat { READS_CUSTOM, FASTQ, CSFASTA };
enum AlignAlgorithm { CPU_DP, GPU_DP };

typedef struct options OPTIONS;

/**
 * \brief This struct is used to maintain options.
 *
 */
struct options {
  // constructor to assign default values
  options();

  // DEBUG INFO
  bool verbose;
  
  // input information
  /**
   * \brief The format of the input genome file
   */
  GenomeFormat genomeFormat;
  /**
   * \brief The format of the input reads file
   */
  ReadsFormat readsFormat;
  /**
   * \brief The path of the input genome file
   */
  string genomeFile;
  /**
   * The path of the input reads file
   */
  string readsFile;

  // ---------------------------------------------------------
  //                   OUTPUT INFORMATION
  // ---------------------------------------------------------
  /**
   * \brief The path of the output genome file (for enabled
   * output)
   */
  string genomeOutputFile;
  /**
   * \brief The path of the output reads file (for enbled output)
   */
  string readsOutputFile;
  /**
   * \brief The path of the alignment file (for enabled alignment)
   */
  string alignOutputFile;

  // ---------------------------------------------------------
  //                 PREPROCESSING INFORMATION
  // ---------------------------------------------------------
  /**
   * \brief The padding for the sequencese.
   *
   * Padding is used especially when GPU algorithms are used
   * (for efficiency reasons). The padding is used to guarantee
   * that the sequence size is a multiple of padding bytes 
   */
  size_t padding;
  /**
   * \brief The number of copies of the (padded) genome sequence.
   *
   * Replicated genome sequences is used especially when GPU
   * algorithms are perfomed.
   */
  size_t genomeCopies;

  // ---------------------------------------------------------
  //                   ALGORITHM INFORMATION
  // ---------------------------------------------------------

  /**
   * \brief Selects the task to be performed
   */
  int task;

  /**
   * \brief Enables the output translation in the custom
   * format
   */
  bool translate;
  /**
   * \brief Enables alignment algorithms and their output
   */
  bool align;
  /**
   * \brief Selects the align algorithm
   */
  AlignAlgorithm alignAlgorithm;
  /**
   * \brief Length of a k-mer (when k-mers are used)
   */
  size_t kmerSize;

  // ---------------------------------------------------------
  //                   PARSING OPERATIONS
  // ---------------------------------------------------------
  /**
   * \brief Prints usage instructions and terminates the
   * eecution.
   *
   * \param os The stream where the instruction will be printed
   * \param name The name of the executable
   * \param exitCode The exit code of the program
   *
   * \sa printOptions()
   */
  void printUsage(ostream& os, const char* name, int exitCode);
  /**
   * \brief Extract options from the input line arguments
   *
   * \param argc The number of input arguments (as \c argc
   * of \c main)
   * \param argv The input options (as \c argv of \c main )
   */
  void parseInputArgs(int argc, char** argv);

  // ---------------------------------------------------------
  //                   INFORMATION PRINTING
  // ---------------------------------------------------------
  /**
   * \brief Prints the options information
   * 
   * \param os The stream where output will bi printed
   * \sa printUsage()
   */
  void printOptions(ostream& os);

};

#endif
