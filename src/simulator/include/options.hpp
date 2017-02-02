#ifndef OPTIONS_H
#define OPTIONS_H

enum class Task { Test, Offline, Online, Oracle, ScoreEst, AlignScore, EditDist, Undefined };

#include <string>


enum OpMode { Test, Offline, Online, Oracle, ScoreEst, AlignScore, EditDist };

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

  std::string inputReference; // -i
  std::string inputSAM; // -S
  std::string outputDistribution; // -D
  std::string outputCDF; //-C

  int approxLevel; // -A

  size_t floatPrecision; // TODO

  Task task; // --task (will substituted mode in -O)
  OpMode mode;  // -O
  int subTask;  // -B
  int optFlags; // -f

  bool pipeline;
  bool online; // -v 

  bool verbose; // -v

  int k; // -k

  static  Options opts;
}
;
void parseArguments(int argc, char** argv);


#ifdef HAVE_BOOST_PROGRAM_OPTIONS

void parseArgumentsBoost(int argc, char** argv);

#endif

#endif
