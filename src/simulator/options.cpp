#include "common.h"

#include <getopt.h>
#include <cstdlib>
#include <iostream>

void printUsage() {
  std::cout << "USAGE\n";
  std::cout << "\tsimulator -N [ref_len] -m [read_len] -M [read_count]\n\n";
}

void printArguments() {
  std::cout << "\t\t    +++++  ARGUMENTS +++++\n\n";
  std::cout << "Reference length (N)         " << Options::opts.N << std::endl;
  std::cout << "Reads length (m)             " << Options::opts.m << std::endl;
  std::cout << "Reads count (M)              " << Options::opts.M << std::endl;
  std::cout << "Base error probability (pe)  " << Options::opts.pe << std::endl;
}

void parseArguments(int argc, char** argv) {

  // default options
  Options::opts.N = 1000000; // size of the reference;
  Options::opts.m = 50;
  Options::opts.M = 10;
  Options::opts.pe = 0.01;
  Options::opts.online = false;
  Options::opts.verbose = false;
  
  char c;
  while ((c = getopt(argc, argv, "N:m:M:e:ohv")) != -1) {
    switch(c) {
    case 'N':
      Options::opts.N = atoi(optarg);
      break;
    case 'M':
      Options::opts.M = atoi(optarg);
      break;
    case 'm':
      Options::opts.m = atoi(optarg);
      break;
    case 'e':
      Options::opts.pe = atof(optarg);
      break;
    case 'o':
      Options::opts.online = true;
      break;
    case 'v':
      Options::opts.verbose = true;
      break;
    case 'h':
      printUsage();
      exit(0);
    default:
      printUsage();
      exit(1);
    }
  }

  if (Options::opts.verbose) {
    printArguments();
  }
}
