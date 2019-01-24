#include <include/common.hpp>
#include <include/options.hpp>

#include <include/log.hpp>

#include <getopt.h>
#include <iostream>
#include <limits>
#include <map>

#ifdef HAVE_BOOST_PROGRAM_OPTIONS

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#endif

Options Options::opts;

std::map<std::string, Task> text_to_op =
  {
    {"test", Task::Test} ,
    {"offline", Task::Offline},
    {"online",  Task::Online},
    {"oracle", Task::Oracle},
    {"score", Task::ScoreEst},
    {"align", Task::AlignScore},
    {"edit", Task::EditDist}
  };

void printUsage() {
  std::cout << std::endl;
  std::cout << "USAGE\n";
  std::cout << "\tsimulator [<OPTIONS>]\n\n";
  std::cout << " OPTIONS\n\n";
  std::cout << "\t-N <length>   Length of reference sequence\n";
  std::cout << "\t-M <length>   Length of reads  [1,2,...]\n";
  std::cout << "\t-m <count>    Number of reads  [1,2,...]\n";
  std::cout << "\t-e <error>    Error probability for base call [0.0,1.0]\n";
  std::cout << "\t-P <prec>     The precision required for task based quantity [0.0, 1.0]\n";
  std::cout << "\t-c <conf>     The confidence (or z_confidence) required (task based)\n";
  std::cout << "\t-d <oscill>   Oscillation tolerated (task based) \n";
  std::cout << "\t-k <number>   Generic numeric parameter to be used based on task\n";  
  std::cout << "\t-i <path>     Path of a fasta file for the reference\n";
  std::cout << "\t-D <path>     Path of a output file for distribution of scores\n";
  std::cout << "\t-C <path>     Path of a output file for CDF (cumulative distribution) of scores\n";
  std::cout << "\t-A <approx>   Set approximation levels (default no approx)\n";
  std::cout << "\t-O <op_mode>  Set operation mode of the simulator [0,4] (run -h for details)\n";
  std::cout << "\t-B <sub_task> Defines the sub task for a given mode\n";
  std::cout << "\t-f <flags>    Flags code to be activated (operation mode and subtask dependent\n";
  std::cout << "\t-t <nthr>     The number of threads to be used (not suppoerted by all tasks)\n";
  std::cout << "\t-p            Outputs on standard out for pipelining\n";
  std::cout << "\t-h            Shows extended help\n";
  std::cout << "\t-v            Activate verbose mode\n";
  std::cout << std::endl;
}

void
printArguments() {
  std::cout << "\t\t    +++++  ARGUMENTS +++++\n\n";
  std::cout << "Reference length (N)         " << Options::opts.N << std::endl;
  std::cout << "Reads length (m)             " << Options::opts.m << std::endl;
  std::cout << "Reads count (M)              " << Options::opts.M << std::endl;
  std::cout << "Base error probability (pe)  " << Options::opts.pe << std::endl;
  std::cout << "Operation mode (O)           "
	    << static_cast<int>(Options::opts.task) << std::endl;
}

void printOperationModeDescription() {
  std::cout << " OPERATION MODES\n\n";
  std::cout << "\t 0 (Test)        Performs several tests of quantities\n";
  std::cout << "\t 1 (Offline)     Offline simulation: all genome and reads generated before scoring (Defualt)\n";
  std::cout << "\t 2 (Online)      Online simulation: genome and reads are generated 'on-demand'\n";
  std::cout << "\t 3 (Oracle)      Simulation of M-1 independent pairs\n";
  std::cout << "\t 4 (ScoreEst)    Tests the approximation E[score] ~= E[num]/E[den]\n";
  std::cout << "\t 5 (AlignScore)  Given a SAM files and a reference evaluate the score\n";
  std::cout << "\t 6 (EditDist)    Performs several experiment on the edit distance\n";
}

void
setDefualtParams() {
   // default options
  Options::opts.N = 1000000; // size of the reference;
  Options::opts.m = 50;
  Options::opts.M = 10;
  Options::opts.pe = 0.01;
  Options::opts.precision = 0.01;
  Options::opts.confidence = 0.95;
  Options::opts.oscillation = 0.5;
  
  Options::opts.empiricalDistributionStep = 100;

  Options::opts.inputReference = "";
  Options::opts.inputSAM = "";
  Options::opts.outputDistribution = "";
  Options::opts.outputCDF = "";
  Options::opts.approxLevel = -1;

  Options::opts.floatPrecision = std::numeric_limits< double >::max_digits10;
  

  Options::opts.task = Task::Test;
  Options::opts.subTask = 0;
  Options::opts.optFlags = 0;
  Options::opts.online = false;
  Options::opts.pipeline = false;
  Options::opts.verbose = false;

}

// ----------------------------------------------------------------------
//                      CUSTOM PARSING FUNCTIONS
// ----------------------------------------------------------------------

Task intToTask(int tCode) {
  Task task = Task::Undefined;
  switch(tCode) {
  case 0:
    task = Task::Test;
    break;
  case 1:
    task = Task::Offline;
    break;
  case 2:
    task = Task::Online;
    break;
  case 3:
    task = Task::Oracle;
    break;
  case 4:
    task = Task::ScoreEst;
    break;
  case 5:
    task = Task::AlignScore;
    break;
  case 6:
    task = Task::EditDist;
    break;

  default:
    task = Task::Undefined;
  }
  return task;
}

std::string taskToString(Task task) {
  switch(task) {
  case(Task::Test): return "Test";
  case(Task::Offline): return "Offline";
  case(Task::Online): return "Online";
  case(Task::Oracle): return "Oracle";
  case(Task::ScoreEst): return "ScoreEst";
  case(Task::AlignScore): return "AlignScore";
  case(Task::EditDist): return "EditDist";
  case(Task::Undefined): return "Undefined";
  default: return "TBD";
  }
}

Task parse_opmode_string(std::string param) {
  for (size_t i = 0; i < param.size(); ++i) {
    param[i] = std::tolower(param[i]);
  }
  if (text_to_op.count(param)) {
    return text_to_op[param];
  }
  return Task::Undefined;
}

// ----------------------------------------------------------------------
//                        ARG PARSING FUNCTIONS
// ----------------------------------------------------------------------


void parseArguments(int argc, char** argv) {

  setDefualtParams();
  
  char c;
  while ((c = getopt(argc, argv, "N:m:M:e:P:c:d:k:i:S:D:C:A:O:B:f:t:phv")) != -1) {
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
    case 'P':
      Options::opts.precision = atof(optarg);
      break;
    case 'c':
      Options::opts.confidence = atof(optarg);
      break;
    case 'd':
      Options::opts.oscillation = atof(optarg);
      break;
    case 'k':
      Options::opts.k = atoi(optarg);
      break;
    case 'i':
      Options::opts.inputReference = optarg;
      break;
    case 'S':
      Options::opts.inputSAM = optarg;
      break;
    case 'D':
      Options::opts.outputDistribution = optarg;
      break;
    case 'C':
      Options::opts.outputCDF = optarg;
      break;
    case 'A':
      Options::opts.approxLevel = atoi(optarg);
      break;
    case 'O':
      Options::opts.task = static_cast<Task>(atoi(optarg));
      break;
    case 'B':
      Options::opts.subTask = atoi(optarg);
      break;
    case 'f':
      Options::opts.optFlags = atoi(optarg);
      break;
    case 'p':
      Options::opts.pipeline = true;
      break;
    case 'v':
      Options::opts.verbose = true;
      break;
    case 't':
      Options::opts.n_threads = atoi(optarg);
    case 'h':
      printUsage();
      printOperationModeDescription();
      std::cout << "\n\n- " << optarg << " -\n\n";
      exit(0);
    default:
      printUsage();
      exit(1);
    }
  }
}


#ifdef HAVE_BOOST_PROGRAM_OPTIONS

void
parseArgumentsBoost(int argc, char** argv) {
  std::string infoArg = "";

  setDefualtParams();
  po::options_description opts_desc("All options");
  // set options
  opts_desc.add_options()

    ("info", po::value<std::string>(&infoArg), "Gives information")
    
    ("help,h", "Show help message") // -h, --help

    ("length,N", po::value<size_t>(&Options::opts.N), // -N, --length
     "Size of the genome or strings in general")

    ("read-length,m", po::value<size_t>(&Options::opts.m), // -m ,--read-length
     "Length of the read")

    ("read-count,M", po::value<size_t>(&Options::opts.M), // -M, --read-count
     "Number of reads")
    
    ("error-prob,e", po::value<double>(&Options::opts.pe), // -e, --error-prob
     "Error probability")

    ("max-error,P", po::value<double>(&Options::opts.precision), // -P, --max-error
     "The maximum error tollerated for estimation (epsilon)")

    ("confidence,c", po::value<double>(&Options::opts.confidence), // - , --confidence
     "Confidence interval as a multiple of sigmas of a normal distribution (delta)")

    ("iterations,k", po::value<int>(&Options::opts.k), // -k, --iterations
     "Number of iterations or samples")

    ("input-reference,i",po::value<std::string>(&Options::opts.inputReference), // -i, --input-reference
     "Reference DNA/RNA input file (fasta format)")

    ("input-sam,S", po::value<std::string>(&Options::opts.inputSAM), // -S, --input-sam
     "Sam alignment input file")

    ("output-density,D", po::value<std::string>(&Options::opts.outputDistribution), // -D, --output-density
     "File were density will be written. If left unspecified, no output will be produced")

    ("output-distribution,C", po::value<std::string>(&Options::opts.outputCDF), // -C, --output-distribution
     "File were distribution will be written. If left unspecified, no output will be produced")

    ("approx-level,A", po::value<int>(&Options::opts.approxLevel), // -A, --approx-level
     "Set the approximation level")

    ("task,O", po::value<std::string>(), "Selects the task to be performed")
    
    ("sub-task,B", po::value<int>(&Options::opts.subTask), // -B, --sub-tast
     "Defines the subtask for the operation mode selected") 

    ("flags,f", po::value<int>(&Options::opts.optFlags), // -f, --flads
     "Sets flags to define behavior of mode / subtask")

    ("pipeline,p", "If set will run in pipeline") // -p, --pipeline

    ("verbose,v", "If set will run in verbose output") // -v, --verbose

    ("threads,t", po::value<size_t>(&Options::opts.n_threads),
     "Specify the number of threads to be used")
    
    ; 
  
  // invoke boost command line parser
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, opts_desc), vm);
  po::notify(vm);

  // When info is requested every other option is ignored
  if (vm.count("info")) {
    std::cout << "Information on " << vm["info"].as<std::string>() << std::endl;
    exit(0);
  }

  // set parameters based on passed arguments
  if (vm.count("help")) {
    std::cout << opts_desc << std::endl;
    printOperationModeDescription();
    std::cout << std::endl;
    exit(0);
  }

  if (vm.count("task")) {
    std::string arg = vm["task"].as<std::string>();
    try {
      Options::opts.task= static_cast<Task>(std::stoi(arg));
    } catch(std::invalid_argument&) {
      logInfo("OpMode String");
      Options::opts.task =
	parse_opmode_string(vm["task"].as<std::string>());
    }
  }

  if (vm.count("pipeline")) {
    Options::opts.pipeline = true;
  }

  if (vm.count("verbose")) {
    Options::opts.verbose = true;
  }

  // DEBUG
  logInfo("Task: " + taskToString(Options::opts.task));


}

#endif


