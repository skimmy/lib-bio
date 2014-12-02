#ifndef _UTIL_IO_
#define _UTIL_IO_

#include<string>

/**
   \file io.hpp
   \brief Collection of utiloty I/O functions
 */

/**
   \fn getFileLength(const std::string& filePath);
   \brief Returns the length in byte of the file with path `filePath`
 */
size_t getFileLength(const std::string& filePath);

template <class T>
void printElements(const T& coll, std::ostream& os, std::string sep = " ", std::string trail = "\n") {
  typename T::const_iterator it;
  for (it = coll.begin(); it != coll.end(); ++it) {
    os << *it << sep;
  }
  os << trail;
}

#endif
