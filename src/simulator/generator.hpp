#ifndef GENERATOR_H
#define GENERATOR_H

#include <string>
#include <queue>

namespace lbio { namespace sim { namespace generator {

struct Read {
  size_t j;  
  std::string r;
  Read(const std::string& s, size_t p) : j(p), r(s) {}
  // !!! WARNING '>' is used to make the priority queue work in
  // !!! ascending order
  bool operator < (const Read& other) const { return this->j > other.j; }
};

void
initRandomGenerator();

void
simpleIIDErrors(std::string& s, double pe);

void
generateIIDGenome(size_t N, char* S);
void
generateIIDString(std::string& s);
void
generateConstantGenome(size_t N, char* S, char b);

void
simulateReadAt(size_t j, size_t m, const char* S, char* r);
void
generateOfflineReads(const std::string& s, std::priority_queue<Read>& reads,
		     size_t m, size_t M, double pe);
Read
generateOnlineRead(char* S, size_t j, size_t m, double pe);
Read
randomRead(size_t m);

void
complementBases(char* S, size_t n);
char
randomMutation(char c);
char
baseComplement(char b);

class IidPairGenerator {
public:

  typedef typename std::pair< std::string, std::string > StringPair;
  
  IidPairGenerator(size_t n_, size_t m_) :
    n {n_}, m {m_}, s1(n_, 'N'), s2(m, 'N') { }
  StringPair operator()() {
    generateIIDString(s1);
    generateIIDString(s2);
    return StringPair {s1, s2};    
  }
  StringPair last_pair() {
    return StringPair {s1, s2};
  }
private:
  size_t n;
  size_t m;
  std::string s1;
  std::string s2;
};
      
} } } // namespaces

#endif
