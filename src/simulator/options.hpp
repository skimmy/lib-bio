#ifndef OPTIONS_H
#define OPTIONS_H

#define HAVE_BOOST


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

  OpMode mode;  // -O
  int subTask;  // -B
  int optFlags; // -f

  bool pipeline;
  bool online; // -v 

  bool verbose; // -v

  int k; // -k

  static  Options opts;
};

void parseArguments(int argc, char** argv);


#ifdef HAVE_BOOST

void parseArgumentsBoost(int argc, char** argv);

#endif

#endif
