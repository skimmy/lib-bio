#ifndef SIM_LOG_H
#define SIM_LOG_H

#include <iostream>

const std::string TermColHeader = "\033[1;95m";
const std::string TermColOkblue = "\033[1;94m";
const std::string TermColOkgreen = "\033[1;92m";
const std::string TermColWarning = "\033[1;93m";
const std::string TermColFail = "\033[1;91m";

const std::string TermColBold = "\033[1m";
const std::string TermColUnderline = "\033[4m";
const std::string TermColEndColor = "\033[0m";


std::string colorifyString(const std::string& plain, const std::string& color);
  

// If debug is disabled or log operations are mapped to no operation
#if defined NDEBUG
#define logError(x)  LBIO_NOOP
#define logWarning(x) LBIO_NOOP
#define logInfo(x) LBIO_NOOP
#else
#define logError(x) _logError(x)
#define logWarning(x) _logWarning(x)
#define logInfo(x) _logInfo(x)
#endif

void setLogStream(std::ostream* stream);

void _logError(const std::string& msg);
void _logWarning(const std::string& msg);
void _logInfo(const std::string& msg);

void fatal_error(const std::string &msg, int exit_code = 1);

#endif
