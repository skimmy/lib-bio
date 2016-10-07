#ifndef ONLINE_H
#define ONLINE_H

// used to represent currenct 'under sequence' segment of a large genome (used
// when simulator is operating in 'online' mode)
struct GenomeSegment {
  char* genome;

  size_t length;
  size_t current_start;
  size_t total_length;
  size_t read_length;
  
  
  GenomeSegment(size_t N, size_t m, size_t l);
  ~GenomeSegment();
};

const size_t MAX_GENOME_SEGMENT_LENGTH = 1 << 20;
void generateFirstGenomeSegment(GenomeSegment& g);
void generateNewGenomeSegment(GenomeSegment& g, size_t keep_pref);

#endif
