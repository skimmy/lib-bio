#include "io.hpp"

#include <sys/stat.h>

size_t getFileLength(const std::string& filePath) {
  size_t n = 0;
  struct stat buffer;
  int status = stat(filePath.c_str(), &buffer);
  if (status == 0) {
    n = (size_t) buffer.st_size;
  }
  return n;
}
