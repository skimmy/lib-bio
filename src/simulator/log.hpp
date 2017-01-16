#ifndef SIM_LOG_H
#define SIM_LOG_H

#include <iostream>

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
