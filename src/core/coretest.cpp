#include "Read.hpp"

#include <iostream>

int main(int argc, char** argv) {
  Read r;
  r.setBases("ACG");
  std::cout << r.getBases() << std::endl;
  return 0;
}
