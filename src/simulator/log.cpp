#include "log.hpp"

#include <iostream>

std::ostream* logStream = &std::cerr;

void
setLogStream(std::ostream* stream) {
  logStream = stream;
}

void
_logError(const std::string& msg) {
  (*logStream) << "\033[1;31mError   \033[0m" << msg << std::endl;
}

void
_logWarning(const std::string& msg) {
  (*logStream) << "\033[1;33mWarning \033[0m" << msg << std::endl;
}

void
_logInfo(const std::string& msg) {
  (*logStream) << "\033[1;34mInfo    \033[0m" << msg << std::endl;
}

void
fatal_error(const std::string &msg, int exit_code)
{
  logError(msg);
  exit(exit_code);
}
