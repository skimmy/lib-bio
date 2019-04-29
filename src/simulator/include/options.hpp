#ifndef OPTIONS_H
#define OPTIONS_H

enum class Task { Test, Offline, Online, Oracle, ScoreEst, AlignScore, EditDist, Undefined };

#include <string>


// options parsing
struct Options {
  size_t N; // -M
  size_t m; // -m
  size_t M; // -M
  double pe; // -e

  size_t empiricalDistributionStep;
  double precision;   // -P
  double confidence;  // -c
  double oscillation; // -d

  std::string alphabet; // -a

  std::string inputReference; // -i
  std::string inputSAM; // -S
  std::string outputDistribution; // -D
  std::string outputCDF; //-C

  int approxLevel; // -A

  size_t floatPrecision; // TODO
  unsigned int seed; // -s

  Task task; // --task (will substituted mode in -O)
  int subTask;  // -B
  int optFlags; // -f
  int weight_scheme; // -W

  bool pipeline;
  bool online; // -v 

  int verbose; // -v
  std::string verboseOutput; // -V

  int k; // -k

  size_t n_threads; // -t, --threads

  static  Options opts;
}
;
void parseArguments(int argc, char** argv);


#ifdef HAVE_BOOST_PROGRAM_OPTIONS

void parseArgumentsBoost(int argc, char** argv);

#endif

#endif
