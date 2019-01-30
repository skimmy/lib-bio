#ifndef GENERATOR_H
#define GENERATOR_H

#include <cmath>
#include <string>
#include <queue>

#include <include/util.hpp>

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
initRandomGenerator(lbio_size_t sigma);

void
simpleIIDErrors(std::string& s, double pe);

void
generateIIDGenome(size_t N, char* S);
void
generateIIDString(std::string& s, std::string alphabet);
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
  
  IidPairGenerator(size_t n_, size_t m_, std::string a) :
    n {n_}, m {m_}, s1(n_, 'N'), s2(m, 'N'), alphabet {a} { }
  StringPair operator()() {
    generateIIDString(s1, alphabet);
    generateIIDString(s2, alphabet);
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
  std::string alphabet;
};

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                           ALPHABET ITERATOR                              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class AlphabetIterator {
private:
  using numeric_type = uint64_t;
  using string_type = std::string;
  

  lbio_size_t size;
  string_type alphabet;
  numeric_type max_value;
  numeric_type current;  
  bool end_iter;  
  
public:
  AlphabetIterator(lbio_size_t n, std::string al)
    : size {n}, alphabet {al},
      max_value { static_cast<numeric_type>(std::pow(al.size(), n)) },
      current {0}, end_iter{false} { }

  AlphabetIterator& operator++() {
    end_iter = (++current >= max_value);
    return *this;
  }

  string_type operator*() {
    string_type out {""};
    numeric_type current_copy {current};
    for (lbio_size_t i = 0; i < size; ++i) {
      out = alphabet[current_copy & 0x3] + out;
      current_copy >>= 2;
    }
    return out;
  }

  bool operator==(const AlphabetIterator& other) {
    return ( end_iter && other.end_iter )
      || ( !end_iter && !other.end_iter && current == other.current);
      
  }

  bool operator!=(const AlphabetIterator& other) {
    return !(*this == other);
  }

  AlphabetIterator end() const {
    static AlphabetIterator _end(size, alphabet);
    _end.end_iter = true;
    return _end;
  }
  
};

class SigmaNIterator {
private:
  std::vector<uint8_t> x;
  std::string alphabet;
  bool end;
  
public:
  SigmaNIterator(lbio_size_t n, std::string Sigma) :
    x(n, 0), alphabet {Sigma}, end {false}
    { }

  SigmaNIterator() : end {true} { }

  SigmaNIterator& operator++() {
    if (!end) {
      int j = x.size()-1;
      while(j>=0 and x[j] == (alphabet.size()-1)) {
	x[j] = 0;
	j--;
      }
      if (j>=0) {
	x[j]++;
      }
      else {
	end = true;	
      }
    }
    return *this;
  }

  std::string operator*() {
    std::string s {""};
    for (auto it = x.begin(); it != x.end(); ++it) {
      s += alphabet[*it];
    }
    return s;
  }

  bool operator==(const SigmaNIterator& other) {
    return (end == other.end);
  }

  bool operator!=(const SigmaNIterator& other) {
    return !(*this == other);
  }
  
};

std::vector<std::string>
all_strings(lbio_size_t n, std::string alphabet);

std::vector<std::pair<std::string, lbio_size_t>>
all_invariant_strings_for_partition(Partition& p, std::string alphabet);

std::vector<std::pair<std::string, lbio_size_t>>
permutation_invariant_strings_with_multiplicity(lbio_size_t n, std::string alphabet);
      
} } } // namespaces

#endif
