#include "sequence.h"
#include "io.h"
#include "util.h"
#include "alignment.h"
#include "quality.h"
#include "adt.h"
#include "tasks.hpp"

#include <vector>
#include <iostream>
#include <fstream>

#include <math.h>

#define DEBUG 1

using namespace std;

#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
using namespace boost;

typedef DNACompressedSymbol DnaSymb;

void test() {
  cout << "-- CPU SMITH WATERMAN TEST --\n";
  string a = "AATGTTA";
  string aa = "ACGGT";
  string aaa = "ACGT";
  string b = "AATGTGACGTTTG";
  //  SmithWatermanDP swa(a.c_str(), a.length(), b.c_str(), b.length());
  SmithWatermanDP swa(a, b);
  swa.enableBacktrack();
  swa.computeMatrix();
  swa.printMatrix();
  MatrixPoint2D maxPos = swa.getGlobalBest();
  cout << "Max Position (" << maxPos.i << ", " << maxPos.j << "): " << swa.getScoreAt(maxPos)  << " " << endl;
  swa.printBacktrackMatrix();

  vector<Read> v;
  Read r;
  r.setBases(a);
  v.push_back(r);
  r.setBases(aa);
  v.push_back(r);
  r.setBases(aaa);
  v.push_back(r);
  vector<Position<int>> alignsVector = alignReads(v,b);
  for(Position<int> p : alignsVector) {
    cout << "(" << p.getSequenceId() << ", " << p.getPosition() << ")" << endl;
  }
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

int gpu_main(int argc, char** argv)
{
  // parse command line arguments
  OPTIONS opts;
  opts.parseInputArgs(argc, argv);
  if(DEBUG) {
    opts.printOptions(cout);
  }

  // load reads and genome
  CompressedReadSet* readSet = loadReads(opts);
  CompressedSequence* compGenome = loadGenome(opts);
  // 'genome' is a pointer to the replicated genome with padding
  // 'reads' is a pointer to the reads aligned as we need
  // const void* genome = compGenome->getSequence();
  // const void* reads = readSet->getSequence();
  // --------------------------
  // --- HERE THE CUDA CODE ---
  // --------------------------

  // --- TEST CODE ---
  if (opts.alignAlgorithm == CPU_DP) {
    test();
  }
  // -----------------

  // save data
  if (opts.translate) {
    readSet->writeToFile(opts.readsOutputFile);
    compGenome->writeToFile(opts.genomeOutputFile);
  }

  // free memory
  if (readSet) {
    delete readSet;
  }
  if (compGenome) {
    delete compGenome;
  }
  
  // vector<Read> reads;
  // FastqFormat fastq;
  // fastq.openFile("/home/skimmy/data/ecoli.solid/ecoli_sample_read1.fastq");

  // while(fastq.hasNextRead()) {
  //   reads.push_back(fastq.getNextRead().trim(48));
  // }
  // vector<Read> otherVector;
  // for (size_t i=1; i <reads.size(); ++i ) {
  //   otherVector.push_back(reads[i]);
  // }

  // size_t i = 0;
  // //  size_t N = 1000000;
  // while( !readsFile.eof() ) {
  //   readsFile >> read;
  //   if (read.length() > 0 ) {
  //     read.trim(48);
  //     reads.push_back(read);
  //     //    cout << read;
  //     i++;
  //   }
  // }
  
  // CompressedReadSet crs(reads);
  // crs.writeToFile("/home/skimmy/data/ecoli.solid/sample_gpu.dat");
  return 0;
}

void byesianInference(double* prior, double* evidence, double* posterior, size_t n) {
  double dotProduct = 0;
  for (size_t i = 0; i < n; ++i) {
    dotProduct += prior[i] * evidence[i];
  }
  for (size_t i = 0; i < n; i++ ) {
    posterior[i] = ( prior[i] * evidence[i] ) / dotProduct;
  }
}

int test_main(int argc, char** argv) {
  OPTIONS opts;
  opts.parseInputArgs(argc, argv);
  vector<CSFastRead> reads;
  vector<CSFastRead> otherVector;
  CSFastFormat cs_in(argv[1],argv[2]);
  while(cs_in.hasNextRead()) {
    otherVector.push_back(cs_in.getNextRead());
  }

  ifstream colors_is(argv[1]);
  ifstream qual_is(argv[2]);
  CSFastRead tmp;
  while( !colors_is.eof() || !qual_is.eof() ) {
    tmp.loadBasesAndQualitiesFromFiles(colors_is, qual_is);
    reads.push_back(tmp);
  }

  FastFormat ff;  
  ff.loadFromFile(argv[3]);
  Reference ref = (Reference)ff;
  // ref.toColors('T');
  
  // size_t k = 10;
  // list< FullyQualifiedSequence< ColorAlphabet > > qualsKmers;
  // //  BOOST_FOREACH( CSFastRead r, reads ) {
  // for ( size_t l = 0; l < 50; ++l) {
  //   list< FullyQualifiedSequence< ColorAlphabet > > tmp = reads[l].getFullyQualifiedKMerList(k);
  //   qualsKmers.insert(qualsKmers.end(), tmp.begin(), tmp.end());
  // }


  // BOOST_FOREACH( FullyQualifiedSequence< ColorAlphabet > fqs, qualsKmers ) {
  //   char * seq = (char*)(fqs.getSequence()->getSequence());
  //   for (size_t j = 0; j < k; ++j){
  //     cout << seq[j] << " (" << (*fqs.getQuality())[seq[j]][j] << ")\t";
  //   }
  //   cout << endl;
  //   //    cout << (char*)(fqs.getSequence()->getSequence()) << (*fqs.getQuality())[] << "\n";
  // }
  // cout << endl;

  size_t K = opts.kmerSize;
  list<KMer> refKMers = ref.getKMerList(K);
  unordered_map<KMer, list< size_t >, sequence_hash, sequence_equality> kmersMap;
  // size_t N = refKMers.size();
  // for (size_t i = 1; i <= N ; ++i) {
  //   kmersMap[refKMers.front()]=-i;
  //   refKMers.pop_front();
  // }
  size_t R = 30000; //reads.size();
  // iterate over all reads and all kmers for each read
  
  for (size_t r = 0; r < R; ++r) {
    list<KMer> readKmers = reads[r].getKMerList(K);
    size_t j = 0;
    BOOST_FOREACH(KMer kmer, readKmers) {
      kmersMap[kmer].push_back( r*j );
      ++j;
    }
  }
  // list<size_t>* maxList = 0;
  // KMer* maxKmer = 0;
  // size_t maxCount = 0;
  // typedef unordered_map<KMer, list< size_t >, sequence_hash, sequence_equality> myMap;
  // for (myMap::iterator iter = kmersMap.begin(); iter != kmersMap.end(); ++iter) {
  //   if ( iter->second.size() > maxCount ) {
  //     maxCount = iter->second.size();
  //     maxList = &(iter->second);
  //     maxKmer = &(iter->first);
  //   }
  // }
  // list< size_t > maxFreqKmerList;
  // typedef pair<KMer, list< size_t > > MyPair;
  // BOOST_FOREACH( MyPair p, kmersMap) {

  //   if (maxFreqKmerList.size() < p.second.size() ) {
  //     maxFreqKmerList = p.second;
  //   }
  // }
  
  // cout << maxList->size() << endl;


 
  return 0;
}

int runTask(int argc, char** argv) {
  // parse arguments
  OPTIONS opts;
  opts.parseInputArgs(argc, argv);
  int task = 1; 
  switch(task) {
  case 1:
    {
      string ref = opts.genomeFile; 
      string reads = opts.readsFile; 
      alignFastqReadsSimpleSW(reads, ref, std::cout, 2, 8);      
      break;
    }
  default:
    return 1;
  }
  return 0;
}

int main(int argc, char** argv) {
  //  gpu_main(argc, argv);
  //test_main(argc,argv);
  //test();
  runTask(argc, argv);
  return 0;
}
