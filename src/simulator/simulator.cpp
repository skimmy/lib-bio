/**
 * This is a simulator developed ad-hoc to allow future optimization regardless
 * the changes to other components of the library (which are designed ti be part
 * of a library rather then efficient stand alon tools.
 */

#include <cstdlib>
#include <ctime>

#include <iostream>

#include "../core/Read.hpp"



char bases[] = {'A', 'C', 'G', 'T'};

void printString(char* s, size_t n) {
  for (int i = 0; i < n; ++i) {
    std::cout << s[i];
  }
}

void generateIIDGenome(size_t N, char* S) {
  for (size_t i = 0; i < N; ++i) {
    S[i] = bases[rand() & 0x3];
  }
}

void simulateReadAt(size_t j, size_t m, char* S, char* r) {
  for (size_t l = 0; l < m; ++l) {
    r[l] = S[j+l];
  }
}

int main(int argc, char** argv) {
  char* ref = NULL;
  char* read = NULL;
  srand(time(NULL));
  
  size_t N = 1000000; // size of the reference
  size_t m = 50; // length of reads
  std::cout << "Reference size: " << N << std::endl;
  std::cout << std::endl;
  std::cout << "\t\t+++++  Starting simulation +++++ " << std::endl;
  
  std::cout << "* Reference generation... ";
  ref = new char[N];
  generateIIDGenome(N,ref);
  std::cout << "[OK]" << std::endl;

  std::cout << "* Read generation...";
  read = new char[m];
  simulateReadAt(rand() % N, m, ref, read);  
  std::cout << "[OK]" << std::endl;
  
  printString(read,m); std::cout << std::endl;
  
  std::cout << "* Cleaning...";
  delete[] read;
  delete[] ref;
  std::cout << "[OK]" << std::endl;

  Read rr;
  return 0;
}
