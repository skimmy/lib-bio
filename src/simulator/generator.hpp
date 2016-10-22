#ifndef GENERATOR_H
#define GENERATOR_H

struct Read {
  std::string r;
  size_t j;
  Read(const std::string& s, size_t p) : j(p), r(s) {}
  // !!! WARNING '>' is used to make the priority queue work in ascending order
  bool operator < (const Read& other) const { return this->j > other.j; }
};

void simpleIIDErrors(std::string& s, double pe);

void generateIIDGenome(size_t N, char* S);
void generateConstantGenome(size_t N, char* S, char b);

void simulateReadAt(size_t j, size_t m, const char* S, char* r);
void generateOfflineReads(const std::string& s, std::priority_queue<Read>& reads);
Read generateOnlineRead(char* S, size_t j);

void complementBases(char* S, size_t n);
char randomMutation(char c);
char baseComplement(char b);



#endif