#include "common.hpp"

#include <iostream>

std::ostream* logStream = &std::cerr;

void
setLogStream(std::ostream* stream) {
  logStream = stream;
}

void
logError(const std::string& msg) {
  (*logStream) << "\033[1;31mError   \033[0m" << msg << std::endl;
}

void
logWarning(const std::string& msg) {
  (*logStream) << "\033[1;33mWarning \033[0m" << msg << std::endl;
}

void
logInfo(const std::string& msg) {
  (*logStream) << "\033[1;34mInfo    \033[0m" << msg << std::endl;
}

