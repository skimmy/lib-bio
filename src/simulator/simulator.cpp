/**
 * This is a simulator developed ad-hoc to allow future optimization regardless
 * the changes to other components of the library (which are designed ti be part
 * of a library rather then efficient stand alon tools.
 */

#include <cstdlib>
#include <ctime>

#include <iostream>
#include <list>


char bases[] = {'A', 'C', 'G', 'T'};
char revBases[128];

void initSimulator() {
  revBases['A'] = revBases['a'] = 0;
  revBases['C'] = revBases['c'] = 1;
  revBases['G'] = revBases['g'] = 2;
  revBases['T'] = revBases['t'] = 3;
  srand(time(NULL));
}

char randomMutation(char c) {
  // int C = revBases[c]; std::cout << C << " ";
  // int r = ((rand() % 3) + 1); std::cout << r << " ";
  
  return bases[(revBases[c] + ((rand() % 3) + 1) ) & 0x3];
}

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

void simpleIIDErrors(std::string& s, double pe) {
  size_t m = s.length();
  for (size_t i = 0; i < m; ++i) {
    double x = (double)rand() / RAND_MAX;
    if (x <= pe) {
      s[i] = randomMutation(s[i]);
    }
  }
}

size_t hammingDistance(const char* s1, const char* s2, size_t m) {
  size_t d = 0;
  for (int i = 0; i < m; ++i) {
    if (s1[i] != s2[i]) {
      ++d;
    }
  }
  return d;
}

size_t hammingDistance(const std::string& s1, const std::string& s2, size_t m) {
  return hammingDistance(s1.c_str(), s2.c_str(), m);
}


int main(int argc, char** argv) {
  char* ref = NULL;
  char* read = NULL;
  initSimulator();

  
  size_t N = 1000000; // size of the reference
  size_t m = 50; // length of reads
  size_t M = 10;
  size_t Nbar = N - m + 1;

  double pe = 0.1;
  
  std::cout << "Reference size: " << N << std::endl;
  std::cout << std::endl;
  std::cout << "\t\t+++++  Starting simulation +++++ " << std::endl;
  
  std::cout << "* Reference generation... ";
  ref = new char[N];
  generateIIDGenome(N,ref);
  std::cout << "[OK]" << std::endl;

  std::cout << "* Read generation...";
  std::list<std::string> reads;
  read = new char[m];
  for (size_t h = 0; h < M; ++h) {
    std::string r(ref + (rand() % Nbar),m);
    std::string rr = r;
    std::cout << r << std::endl;
    simpleIIDErrors(r,pe);
    std::cout << r << std::endl << hammingDistance(r,rr,m) << std::endl  << std::endl;
    reads.push_front(r);
  }
  std::cout << "[OK]" << std::endl;
  
  
  std::cout << "* Cleaning...";


  delete[] read;
  delete[] ref;
  std::cout << "[OK]" << std::endl;

  return 0;
}
