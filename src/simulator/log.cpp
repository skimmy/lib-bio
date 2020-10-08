#include <include/log.hpp>

#include <iostream>

std::ostream *logStream = &std::cerr;

void
setLogStream(std::ostream *stream) {
    logStream = stream;
}

void
_logError(const std::string &msg) {
    (*logStream) << "\033[1;31mError   \033[0m" << msg << std::endl;
}

void
_logWarning(const std::string &msg) {
    (*logStream) << "\033[1;33mWarning \033[0m" << msg << std::endl;
}

void
_logInfo(const std::string &msg) {
    (*logStream) << "\033[1;34mInfo    \033[0m" << msg << std::endl;
}

void
_logDebug(const std::string &msg) {
    (*logStream) << TermColDebug + "Debug  \033[0m" << msg << std::endl;
}

void
fatal_error(const std::string &msg, int exit_code) {
    logError(msg);
    exit(exit_code);
}


namespace lbio {
namespace sim {
namespace log {

std::string
colorifyString(const std::string &plain, const std::string &color) {
    return (color + plain + TermColEndColor);
}

std::string
make_bold(const std::string &plain) {
    return (TermColBold + plain + TermColEndColor);
}

std::string
debug_string(std::string head, std::string msg) {
    return (TermColBold + head + TermColEndColor + msg);
}

}
}
} // namespaces
