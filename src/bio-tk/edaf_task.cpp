#include "edaf_task.hpp"

#include <util/str_util.hpp>
#include <../io.h>

#include <algorithms/hamming_distance.hpp>
#include <algorithms/edit_distance.hpp>
#include <algorithms/d2_distance.hpp>

#include <fstream>
#include <iterator>

void
edaf_task(std::map<std::string, std::string> params) {
  lbio_size_t L_max = std::atoi(params["L"].c_str());
  lbio_size_t k = std::atof(params["k"].c_str());
  int verbosity = std::atof(params["verbosity"].c_str());  
  std::string files_opts = params["fastq"];
  std::vector<std::string> files =
    lbio::split(files_opts.cbegin(), files_opts.cend(), ',');
 
  std::ifstream file_1(files[0]);
  std::ifstream file_2(files[1]);

  std::istream_iterator<FastqRead> it_1(file_1);
  std::istream_iterator<FastqRead> it_2(file_2);
  std::istream_iterator<FastqRead> it_end;
  
  lbio::edit_distance_wf<std::string> edit(L_max,L_max);
   while ( (it_1 != it_end) and (it_2 != it_end) ) {
     FastqRead r1 {*it_1};
     FastqRead r2 {*it_2};
     std::string s1 = r1.getBases();
     std::string s2 = r2.getBases();
     lbio_size_t d = lbio::hamming_distance(s1.cbegin(), s1.cend(), s2.cbegin());
     lbio_size_t e = edit.compute(s1, s1.size(), s2, s2.size());
     double D2 = lbio::D2_star<std::string>(k, s1.cbegin(), s1.cend(),
					    s2.cbegin(), s2.cend(),
					    [] (std::string c) { return 0.25; });
     if (verbosity > 0) {
       std::cout << s1 << "\n" << s2 << "\n";
     }
     std::cout<< d << "\t" << e << "\t" << D2 << "\n";
     ++it_1;
     ++it_2;
   }
}
