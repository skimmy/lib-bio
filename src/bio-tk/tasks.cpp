#include <cmath>
#include <vector>
#include <fstream>
#include <iterator>
#include <thread>
#include <unordered_map>
#include <ctime>
#include <algorithm>
#include <random>

#include <util/io_helper.hpp>
// 2D matrix
#include <structures/matrix.hpp>
// prob utilities
#include <prob/generator.hpp>
#include <prob/probability.hpp>


#include "generate_task.hpp"

#include "tasks.hpp"
#include "../core.h"
#include "../io.h"
#include "../algorithms.h"


double rand_double() {
  
  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_real_distribution<double> _dist(0.0,1.0);
  return _dist(gen);
}

//#include "util/io.hpp"

// This utilitu function returns the bytes remaining until the end of the file
// is reached. It should be probably made publicly availale in a separate util
// section of the library (e.g. util/io.h)
size_t bytesRemainingToTheEnd(std::ifstream& ifs) {
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
    std::cout << seq::NumericKMer(p.first, k) << " " << p.second << std::endl;
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
      seq::NumericKMer nkmer(kmer);
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
  // if 0 or other wrong values are set for T, then set it to 1 (i.e. single thread)
  T = (T > 1) ? T : 1;

  // define output (either a file or the standard out)
  std::ofstream outFileStream;
  if (!out.empty()) {
    outFileStream.open(out, std::ios::out);
  }
  std::ostream& outStream = (out.empty()) ? std::cout : outFileStream;
  time_t beginTime, endTime;

  // "split" input file into T threads (if multithreaded is required)
  if (T > 1) {
  }
  
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
  std::cout << " Done (" << M << " reads in " << elapsed  << " sec, "
	    << rate << " reads/sec)"  << std::endl;
}

void
task_read_statistics(const std::string& reads, const string& w_dir,
		     const string& prefix) {

  using BaseQualPair = std::pair<char,int>;
  using IntIntMap = std::map<int,int>;
  IntIntMap lengths;
  std::map<BaseQualPair,lbio_size_t> base_qual_freq;
  
  //IntIntMap qualities;
  //  CharIntMap bases;
  
  lbio_size_t total_pair = 0;
  // for each read in the stream
  std::ifstream ifs(reads, std::ifstream::in);
  for (std::istream_iterator<FastqRead> it(ifs);
       it != std::istream_iterator<FastqRead>(); ++it) {
    FastqRead read {*it};
    lengths[read.length()]++;
    lbio_size_t read_len = read.length();
    std::string bs = read.getBases();
    std::string qs = read.getQualities();
    for (lbio_size_t i = 0; i < read_len; ++i) {
      base_qual_freq[std::make_pair(bs[i], static_cast<int>(qs[i]-33))]++;
      total_pair++;
    }    
  }
  ifs.close();

  // save length stat file
  // [len]\t[count]
  std::string filePath{ w_dir + prefix + "len.btk" };
  std::ofstream lenFile(filePath, std::ofstream::out);
  for (auto key : lengths) {
    lenFile << key.first << "\t" << key.second << "\n";
  }
  lenFile.close();

  // save on file (Base,Qual) frequency
  std::string pair_file_path { w_dir + prefix + "base_qual.btk" };
  std::ofstream pairs_file(pair_file_path, std::ofstream::out);
  std::for_each(base_qual_freq.cbegin(), base_qual_freq.cend(),
		[&pairs_file](std::pair<BaseQualPair, lbio_size_t> pair) {
		  pairs_file << pair.first.first << "\t" << pair.first.second 
			     << "\t" << pair.second << "\n";
		});
  pairs_file.close();
}

void
task_generate(std::map<std::string,std::string> gen_params) {
  using IIDCharSampler = lbio::IIDSampler<lbio::DiscreteProbability<char>>;
  using CharDoubleDistr = lbio::DiscreteProbability<char,double>;
  std::map<char,double> equal_prob_bases =
    {
      { 'A', 0.25},
      { 'C', 0.25},
      { 'G', 0.25},
      { 'T', 0.25}
    };
  CharDoubleDistr _distr(equal_prob_bases.cbegin(), equal_prob_bases.cend());
  IIDCharSampler sampler(_distr);
  
  lbio_size_t _N = atoi(gen_params["N"].c_str());
  lbio_size_t _L = atoi(gen_params["L"].c_str());
  std::ofstream _fastq(gen_params["fastq"]);
  
  // reads from reference
  if (gen_params.count("reference") > 0 && !gen_params["reference"].empty()) {
    // when reference length G is given, reference must be generated
    if (gen_params.count("G") > 0) {
      // generate reference
      lbio_size_t _G = lbio::from_string<lbio_size_t>(gen_params["G"]);
      std::string ref = generate_iid_bases(_G, sampler);
      // save reference to file
      std::ofstream ref_os(gen_params["reference"]);
      // TODO: Divide reference into lines
      ref_os << "> G=" << _G << "\n" << ref << "\n";
      ref_os.close();
      std::uniform_int_distribution<lbio_size_t> position_gen(0,_G-1-_L);
      auto& gen = lbio::global_rand_generator<std::mt19937>();
      // reads from reference
      for (lbio_size_t i = 0; i < _N; ++i) {
	lbio_size_t j = position_gen(gen);
	std::string read = ref.substr(j, _L);

	std::string quals(_L,'!');
	std::stringstream read_stream {};
	read_stream << "> " << i << "pos=" << j << "\n" << read << "\n+\n" << quals << "\n";
    
	FastqRead r;
	read_stream >> r;
	_fastq << r;
      }
    }
    return;
  }
  
 

  
  for (size_t _i = 0; _i < _N; ++_i) {
    std::string seq = generate_iid_bases(_L, sampler);
    // if pattern is in config insert it with specific probability
    if (gen_params.count("motif") and gen_params.count("p_motif")) {
      std::string motif = gen_params["motif"];
      double p_motif = lbio::from_string<double>(gen_params["p_motif"]);
      // implant motif
      for (size_t _j = 0; _j < seq.size() - motif.size(); ++_j) {
	if (rand_double() < p_motif) {
	  seq.replace(_j, motif.size(), motif);
	  _j += motif.size();
	}
      }
    }
    std::string quals(_L,'!');
    std::stringstream read_stream {};
    read_stream << "> " << _i << "\n" << seq << "\n+\n" << quals << "\n";
    
    FastqRead r;
    read_stream >> r;
    _fastq << r;
  }
}
