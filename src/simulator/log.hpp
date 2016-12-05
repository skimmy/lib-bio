#ifndef SIM_LOG_H
#define SIM_LOG_H

#include <iostream>

void
setLogStream(std::ostream* stream);

void
logError(const std::string& msg);

void
logWarning(const std::string& msg);

void
logInfo(const std::string& msg);

#endif
