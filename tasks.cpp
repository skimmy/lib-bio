#include <cmath>
#include <vector>
#include <fstream>
#include <thread>
#include <unordered_map>
#include <ctime>

#include "tasks.hpp"
#include "sequence.h"
#include "io.h"
#include "alignment.h"
#include "algorithms.h"

#include "util/io.hpp"

// This utilitu function returns the bytes remaining until the end of the file
// is reached. It should be probably made publicly availale in a separate util
// section of the library (e.g. util/io.h)
size_t bytesRemainingToTheEnd(ifstream& ifs) {
  std::streampos tmp = ifs.tellg();
  ifs.seekg(0, ifs.end);
  size_t totalLength = ifs.tellg();
  ifs.seekg(tmp);
  return (totalLength - ifs.tellg());
}

void alignSmithWaterman(std::vector<Read>* reads, const Reference* ref, 
			std::vector<ScoredPosition<int,int> >* aligns, int indexOffset) {
  int nReads = reads->size();
  string refBases((char*)ref->getSequence());
  std::cout << "Aligning " << nReads << " reads" << std::endl;
  for (int i = 0; i < nReads; ++i) {    
    SmithWatermanDP sw((*reads)[i].getBases(), refBases);
    sw.computeMatrix();
    MatrixPoint2D maxP = sw.getGlobalBest();
    aligns->push_back(ScoredPosition<int, int>((i + indexOffset),maxP.j, sw.getScoreAt(maxP)));
  }
}



std::vector<ScoredPosition<int, int> > alignFastqReadsSimpleSW(const string& readsPath, const string& referencePath, 
							       std::ostream& output, uint64_t nThreads, size_t nReads) {
  
  std::cout << "-- Smith Waterman alignment --" << std::endl;
  std::cout.flush();
  
  std::vector<ScoredPosition<int,int> > aligns;
						   // some constants
  uint64_t T = (nThreads > 0) ? nThreads : 1;
  uint64_t M = (nReads > 0) ? nReads : 0;

  std::cout << " ************* " << T << "  " <<  M << "**************"  << std::endl;
  std::cout.flush();
	       
  // open files
  std::ifstream readsIn(readsPath, std::ifstream::in);
  //  std::ifstream refIn(referencePath, std::ifstream::in);  

  // list of all threads (use later for joining)  
  std::vector<std::thread> threads;
  std::vector<std::vector<ScoredPosition<int,int> >*> threadAligns;

  // read input reference...
  std::cout << "    Loading reference..." << std::endl;
  std::cout.flush();
  FastFormat fast;
  fast.loadFromFile(referencePath);
  Reference ref = (Reference)fast;  

  // ...and reads
  std::cout << "    Loading reads..." << std::endl;
  std::cout.flush();
  if (nReads >= 0) {
    // CASE 1: A specific total number of reads has been specified or
    uint64_t readsPerThread = (size_t)ceil( ((double)(M)) / ((double)(T)) );
    uint64_t actualReads = 0;    
    // repeat for each thread
    for (uint64_t t = 0; t < T; ++t) {
      std::vector<Read>* reads = new std::vector<Read>();
      actualReads = 0;
      while( (!readsIn.eof()) &&  (actualReads < readsPerThread) ) {	
	FastqRead r;
	readsIn >> r;
	reads->push_back(r);
      	actualReads++;
      }
      // create the aligns vector for the next starting thread
      std::vector<ScoredPosition<int,int> >* alignsVector = new std::vector<ScoredPosition< int,int > >();
      threadAligns.push_back(alignsVector);
      threads.push_back(std::thread(alignSmithWaterman, reads, &ref, alignsVector, t * readsPerThread));
    }


  } else {
    // CASE 2: we need to estimate the amount of reads to be assigned to each thread
    double fraction = 0;
    double D = 1.0 / ((double)T);
    // repeat for each thread
    for (uint64_t t = 0; t < T; ++t) {
      while( (!readsIn.eof()) && (fraction < D)) {
      	fraction += 0.1;
      }
      //      threads.push_back(std::thread(alignSmithWaterman,reads,&ref));
    }
    
  }

  // synchronize all threads
  for (auto& thr : threads) {
    thr.join();
  }
  aligns.clear();
  for (std::vector<ScoredPosition<int,int> >* v : threadAligns) {
    std::cout << "++++ " << v->size() <<  std::endl;
    aligns.insert(aligns.end(), v->begin(), v->end());
    delete v;
  }
  

  output << "*** Alignments: " << std::endl;
  for (ScoredPosition< int,int > p : aligns) {
    output << p.getSequenceId() << "\t" << p.getPosition() << " (" << p.getScore() << ")" << std::endl;
  }

  return aligns;
}

/**************************** K-SPECTRUM FUNCTIONS ****************************/
void taskComputeKSpectrum(size_t k, const string& referenceFile) {
  FastFormat fast;
  fast.loadFromFile(referenceFile);
  Reference ref = (Reference) fast;
  std::unordered_map< uint64_t, uint64_t > index = spectrumAsIntMap(ref, k);
  for (std::pair< uint64_t, uint64_t > p : index) {
    std::cout << NumericKMer(p.first, k) << " " << p.second << std::endl;
  }
}

void taskMapReadsKmers(const string& reference, const string& reads, size_t k, const string& out) {
  std::cout << "-------------------- Reads Mapping --------------------" << std::endl;
  // open files
  std::ofstream outFileStream;
  if (!out.empty()) {
    outFileStream.open(out, std::ios::out);
  }
  std::ostream& outStream = (out.empty()) ? std::cout : outFileStream;
  FastFormat refFast;
  std::cout << "Loading reference (" << reference << ")...";
  std::cout.flush();
  refFast.loadFromFile(reference);
  std::cout << " Ok!" << std::endl;
  Reference ref = refFast.toReference();
  
  std::ifstream readsStream(reads, std::ios::in);

  // compute index for the reference
  std::cout << "Index creation...";
  std::cout.flush();
  NumericKmerIndex index = kmersMapping(ref, k);
  std::cout << " Ok!" << std::endl;
  // scan reads and find mappings
  size_t read_index = 0;
  size_t kmer_index = 0;
  std::cout << "Reads mapping...";
  std::cout.flush();
  while(!readsStream.eof()) {
    FastqRead r;
    readsStream >> r;
    list< KMer > kmers = r.getKMerList(k);
    kmer_index = 0;
    for (KMer kmer : kmers) {
      NumericKMer nkmer(kmer);
      std::unordered_map< uint64_t, std::list< size_t > >::const_iterator it = index.find((uint64_t)nkmer);
      outStream << read_index << ":" << kmer_index << " "  << ((it == index.end()) ? (int64_t)(-1) :  (int64_t)it->second.front() ) << std::endl;
      kmer_index++;
    }
    read_index++;
  }
  std::cout << "Ok!" << std::endl << "-------------------- DONE --------------------" << std::endl;
}


// *****************************************************************************
// *                                   KSCORE                                  *
// *****************************************************************************
// Computes for each reads of the input set the kmerscore and the number of errors
// (see algorithms/kmerscore) against the reference sequence
void taskKmerScoreReads(const string& reference, const string& reads, size_t k, const string& out, size_t T) {

  // define output (either a file or the standard out)
  std::ofstream outFileStream;
  if (!out.empty()) {
    outFileStream.open(out, std::ios::out);
  }
  std::ostream& outStream = (out.empty()) ? std::cout : outFileStream;
  time_t beginTime, endTime;
  
  std::cout << "-------------------- Reads Mapping --------------------" << std::endl;
  // open files 
  std::cout << "Loading reference (" << reference << ")...";
  std::cout.flush(); // in case of stuck code we see at which point
  FastFormat refFast(reference);  
  Reference ref = refFast.toReference();
  std::cout << " Done (" << ref.getSequenceLength() << " bases)" << std::endl;

  // open a stream for reading reads file
  std::ifstream readsStream(reads, std::ios::in);

  // compute the index for the reference
  std::cout << "Index creation...";
  std::cout.flush();
  time(&beginTime);
  NumericKmerIndex index = kmersMapping(ref, k);
  time(&endTime);  
  std::cout << " Ok! (" << difftime(endTime, beginTime) << " sec)" << std::endl;
  std::cout << "Calculating scores...";
  std::cout.flush();
  size_t M = 0;
  time(&beginTime);
  while(!readsStream.eof()) {
    FastqRead r;
    readsStream >> r;
    if (r.getSequenceLength() == 0) {
      continue;
    }
    KmersMap map = extractKmersMapPosition(r, index, k);
    std::vector< uint64_t > scoreVector = kmerScoreVector(map, k);
    KmerScoreType score = scoreForVector(scoreVector, k);
    size_t errors = kmerErrorCount(map, k);
    outStream << score << " " << errors << std::endl;
    M++;
  }
  readsStream.close();
  time(&endTime);
  double elapsed = difftime(endTime, beginTime);
  double rate = M / elapsed;
  std::cout << " Done (" << M << " reads in " << elapsed  << " sec, " << rate << " reads/sec)"  << std::endl;
}
