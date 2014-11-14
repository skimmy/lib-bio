#include "sequence.h"
#include "io.h"
#include "util.h"
#include "alignment.h"
#include "quality.h"
#include "adt.h"
#include "tasks.hpp"
#include "algorithms.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>

#include <math.h>

#define DEBUG 1

using namespace std;

#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
using namespace boost;

typedef DNACompressedSymbol DnaSymb;

void test() {
  // -------------------
  // SMITH WATERMAN TEST
  // -------------------
  // cout << "-- CPU SMITH WATERMAN TEST --\n";
  // string a = "AATGTTA";
  // string aa = "ACGGT";
  // string aaa = "ACGT";
  // string b = "AATGTGACGTTTG";
  // //  SmithWatermanDP swa(a.c_str(), a.length(), b.c_str(), b.length());
  // SmithWatermanDP swa(a, b);
  // swa.enableBacktrack();
  // swa.computeMatrix();
  // swa.printMatrix();
  // MatrixPoint2D maxPos = swa.getGlobalBest();
  // cout << "Max Position (" << maxPos.i << ", " << maxPos.j << "): " << swa.getScoreAt(maxPos)  << " " << endl;
  // swa.printBacktrackMatrix();

  // vector<Read> v;
  // Read r;
  // r.setBases(a);
  // v.push_back(r);
  // r.setBases(aa);
  // v.push_back(r);
  // r.setBases(aaa);
  // v.push_back(r);
  // vector<Position<int>> alignsVector = alignReads(v,b);
  // for(Position<int> p : alignsVector) {
  //   cout << "(" << p.getSequenceId() << ", " << p.getPosition() << ")" << endl;
  // }

  // -----------------
  // NUMERIC KMER TEST
  // -----------------
  // Reference ref("TACGGTGGTCTAA");
  // size_t k = 2;
  // NumericKMer kmer(ref);
  // std::cout << ref << " (" << k << ")" << std::endl;
  // std::unordered_map< uint64_t, std::list< size_t> > index = kmersMapping(ref, k);
  // for (pair< uint64_t, std::list< size_t > > p : index) {
  //   std::cout << NumericKMer(p.first, k) << " --> [ ";
  //   for (size_t i : p.second) {
  //     std::cout << i << " ";
  //   }
  //   std::cout << "]" << std::endl;      
  // }

  // taskMapReadsKmers("/home/skimmy/filtering/data/ecoli.fasta",
  // 		    "/home/skimmy/filtering/data/ecoli.sample.custom.fastq",
  // 		    14, "/home/skimmy/Temporary/ecoli_kmers.out");
}

// TODO
// 2. Create minimal code to perform gpu alignment

/**
 * This function create m copies (shifted on the first base by
 * 0, 1, ... , m-1 bases) of n base long sequence. If n is greater
 * then the real size of sequence, the remaining positions are
 * filled with zeros, otherwise sequence is truncated to the
 * first n bases.
 * @param sequence The original sequence
 * @param n the length of the single sequence
 * @param m the number of shifted copies
 * @return A CompressedSequence dynamically allocated and containing
 * the m shifted copies of the n bases long sequence
 */
CompressedSequence* padAndShiftGenome(const string& sequence, size_t n, size_t m) {
  CompressedSequence* cs = new CompressedSequence(n*m,4);
  size_t N = (size_t)min(sequence.length(), n);
  size_t k = 0;
  for (size_t i = 0; i < m; i++) {
    for (size_t j = m; j < N; j++) {
      cs->setElementAt(k,DnaSymb::IupacToNumber(sequence[j]));
      k++;
    }
  }
  return cs;
}

/**
 * This functions returns a CompressedReadSet object initialized 
 * with the set of reads obtained from the file passed as
 * command line argument and stored in the opt structure. The
 * returned reference is NULL if the type of input file can not
 * be determined, otherwise it points to a dynamically allocated
 * object which must be freed by the programmer.
 *
 * @param otps The structure with parsed command line options
 * @return A CompressedReadSet object upon success NULL otherwise
 */
CompressedReadSet* loadReads(const OPTIONS& opts) {
  CompressedReadSet* crs = new CompressedReadSet();
  switch(opts.readsFormat) {
  case READS_CUSTOM:
    crs->loadFromFile(opts.readsFile);
    break;
  case FASTQ:
    break;
  case CSFASTA:
    break;
  default:
    delete crs;
    return 0;
  }
  return crs;
}

CompressedSequence* loadGenome(const OPTIONS& opts) {
  CompressedSequence* sequence = new CompressedSequence();
  switch (opts.genomeFormat) {
  case GENOME_CUSTOM:
    sequence->loadFromFile(opts.genomeFile);
    break;
  case FAST: {
    delete sequence;
    FastFormat ff;
    ff.loadFromFile(opts.genomeFile);
    size_t n = (size_t) (opts.padding * ceil((double)ff.getSequence().length() / (double)opts.padding));
    return padAndShiftGenome(ff.getSequence(),n, opts.genomeCopies);
  }
  default:
    delete sequence;
    return 0;
  }
  return sequence;
}

int runTask(int argc, char** argv) {
  // parse arguments
  OPTIONS opts;
  opts.parseInputArgs(argc, argv);
  #ifndef DEBUG
    opts.printOptions(std::cout);
  #endif
  int task = opts.task;
  switch(task) {
  case 0:
    {      
      std::cout << std::endl << "'nop' task selected\n\tno operation has been performed\n" << std::endl;
      break;
    }
  case 1: // Smith Waterman alignment
    {
      string ref = opts.genomeFile; 
      string reads = opts.readsFile; 
      alignFastqReadsSimpleSW(reads, ref, std::cout, 2, 8);      
      break;
    }
  case 2: // k-spectrum calculation
    {
      size_t k = opts.kmerSize;
      string ref = opts.genomeFile;
      taskComputeKSpectrum(k, ref);
      break;
    }
  case 3:
    {
      // k-mer mapping
      size_t k = opts.kmerSize;
      string ref = opts.genomeFile;
      string reads = opts.readsFile;
      string out = opts.alignOutputFile;
      taskMapReadsKmers(ref, reads, k, out);
      break;
    }
  default:
    std::cout << "Unrecognized operation" << std::endl;
    return 1;
  }
  return 0;
}

int main(int argc, char** argv) {
  runTask(argc, argv);
  //test();
  
  
  return 0;
}
