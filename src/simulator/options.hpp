#ifndef OPTIONS_H
#define OPTIONS_H

enum OpMode { Test, Offline, Online, Oracle, ScoreEst, AlignScore };

// options parsing
struct Options {
  size_t N;
  size_t m;
  size_t M;
  double pe;

  size_t empiricalDistributionStep;

  std::string inputReference;
  std::string inputSAM;
  std::string outputDistribution;
  std::string outputCDF;

  int approxLevel;

  OpMode mode;

  bool pipeline;
  bool online;

  bool verbose;
  bool test;

  static  Options opts;
};

void parseArguments(int argc, char** argv);

#endif
