#include "common.hpp"

#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>

#define TEST_GENOME_LENGTH 65536
#define TEST_READ_LENGTH 100

const int MC_SAMPLES = 2 << 22;

void
testUtils() {
  std::cout << "\nUTILITY FUNCTIONS/CLASSES TESTS" << std::endl;

  std::cout << "* Geometric progression\n";
  GeometricProgression<int> geom(2,1);
  for (size_t i = 0; i < 10; ++i) {
    std::cout << geom.getNext() << std::endl;
  }

  std::cout << "* Linear Progression\n";
  LinearProgression<double> lin(0.5,0);
  for (size_t i = 0; i < 10; ++i) {
    std::cout << lin.getNext() << std::endl;
  }
}

void
testProbFunctions() {
  std::cout << "\nPROBABILITY FUNCTIONS TESTS\n\n";

  std::cout << "* Median tests\n";
  int f_sym [5] = {1, 3, 4, 3, 1};
  size_t median_idx = medianFromFrequency<int>(f_sym,5);
  std::cout << "  Median for {1, 3, 4, 3, 1}\t\t"
	    << f_sym[median_idx] << " i: " << median_idx << "\n";

  double f_skew [7] = {0.5, 0.5, 0.5, 1, 1, 2, 0.5};
  median_idx = medianFromFrequency<double>(f_skew, 7);
  std::cout << "  Median for {0.5, 0.5, 0.5, 1, 1, 2, 0.5}\t\t"
	    <<  f_skew[median_idx]  << " i: " << median_idx << "\n";

  double f_0 [3] = {0.5, 0.25, 0.25};
  median_idx = medianFromFrequency<double>(f_0, 3);
  std::cout << "  Median for {0.5, 0.25, 0.25}\t\t"
	    <<  f_0[median_idx]  << " i: " << median_idx << "\n";  

}

void
testSampleEstimators() {

  std::cout << "\nESTIMATORS TESTS\n\n";

  std::cout << "Simple test\n Samples: [1, 2, -2, 4, 2, 2]\n";
  size_t k = 6;
  int x[] = {1, 2, -2, 4, 2, 2};
  SampleEstimates est = estimatesFromSamples<int>(x,k);
  std::cout << "* Sample Mean:   " << est.sampleMean << std::endl;
  std::cout << "* Sample Var:    " << est.sampleVariance << std::endl;

  std::cout << "\nConstant test\n Samples = [0.5, 0.5, 0.5]\n";
  double y[] = {0.5, 0.5, 0.5};
  k = 3;
  est = estimatesFromSamples<double>(y,k);
  std::cout << "* Sample Mean:   " << est.sampleMean << std::endl;
  std::cout << "* Sample Var:    " << est.sampleVariance << std::endl;

  int z[128];  
  for (int i = 0; i < 128; ++i) {
    z[i] = (rand() % 65) - 32;
  }
  std::cout << "\nRandom test\n 128 samples U[-32,32]\n";
  est = estimatesFromSamples<int>(z,128);
  std::cout << "* Sample Mean:   " << est.sampleMean << std::endl;
  std::cout << "* Sample Var:    " << est.sampleVariance << std::endl;

}

void testApproximatedExpectedScore() {
  for (int s = 0; s < Options::opts.m + 1; ++s) {
    std::cout << s << "\t" << approximatedScore(s) << "\n";
  }
}

// tests the value of p_eq by simulating 
void testPeq() {

  double pe = Options::opts.pe;
  double jointProb[4][4];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j <4; j++) {
      jointProb[i][j] = 0;
    }    
  }
  // joint probability distribution (R1 = x, R2 = y | R0 = tBase)
  double x = 0;
  size_t i,j;
  char R1,R2;
  for (int i = 0; i < MC_SAMPLES; ++i) {
    char tBase = bases[rand() & 0x03];
    R1 = R2 = tBase;
    x = (double)rand() / RAND_MAX;
    if (x < pe) {
      R1 = randomMutation(tBase);
    }
    x = (double)rand() / RAND_MAX;
    if (x < pe) {
      R2 = randomMutation(tBase);
    }
    jointProb[revBases[R1]][revBases[R2]]++;
  }
  std::cout << "\n";
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j <4; j++) {
      std::cout << jointProb[i][j] / (double)MC_SAMPLES << '\t';
    }
    std::cout << '\n';
  }
}


void testLookupTables() {
  clearUtil();
  initUtil();
  std::cout << "LOOKUP TABLES TEST\n";
  std::cout << "4^{-(m-s)}\n";
  for (size_t i = 0; i <= Options::opts.m; ++i) {
    std::cout << "s: " << i << "\t" << power4_lookup[i] << std::endl;
  }
}

void testScoreFunctionNoError() {
  // Test (i) two reads with complete prefix suffix overlap   
  std::cout << "NO ERRORS\n";
  std::cout << "(i)\n";

  Options::opts.pe = 0;
  Options::opts.m = TEST_READ_LENGTH;
  Options::opts.N = TEST_GENOME_LENGTH;
  
  char S[TEST_GENOME_LENGTH];
  generateConstantGenome(TEST_GENOME_LENGTH, S, 'A');
  char r1[TEST_READ_LENGTH];
  char r2[TEST_READ_LENGTH];
  
  simulateReadAt(0,TEST_READ_LENGTH,S,r1);
  simulateReadAt(0,TEST_READ_LENGTH,S,r2);

  for (size_t i = 1; i <= TEST_READ_LENGTH; ++i) {
    std::cout << i << '\t' << score(std::string(r1), std::string(r2), i) << std::endl;
  }
  

  // TODO this test is 'trivial' leave as last one
  // Test (ii) two reads with no prefix suffix overlap
  std::cout << "(ii)\n";

  Options::opts.pe = 0;
  Options::opts.m = TEST_READ_LENGTH;
  Options::opts.N = TEST_GENOME_LENGTH;
  
  simulateReadAt(0,TEST_READ_LENGTH,S,r1);
  simulateReadAt(0,TEST_READ_LENGTH,S,r2);


  // Test (iii) two reads with exactly one prefix suffix overlap
  std::cout << "(iii)\n";

  size_t m = TEST_READ_LENGTH;
  simulateReadAt(0,TEST_READ_LENGTH,S,r1);
  r1[m-1] = baseComplement(r1[m-1]);
  r1[0] = 'G';
  // test all overlaps except exact one
  for (size_t s = 1; s < m; ++s) {
      simulateReadAt(0,TEST_READ_LENGTH,S,r2);
      r2[s-1] = r1[m-1];
      std::cout << s << '\t' << score(std::string(r1), std::string(r2), s) << std::endl;
  }

}

void testScoreFunctionError() {
  std::cout << "ERRORS\n";
  // Test (i) two reads with complete prefix suffix overlap
  std::cout << "(i)\n";

  // Test (ii) two reads with no prefix suffix overlap
  std::cout << "(ii)\n";

  // Test (iii) two reads with exactly one prefix suffix overlap
  std::cout << "(iii)\n";
}

void testScoreFunction() {
  std::cout << "* Score function tests\n";
  std::cout << std::endl;
  testScoreFunctionNoError();
  std::cout << std::endl;
  testScoreFunctionError();

}

void testEditDistance() {
  std::cout << "* Edit distance tests\n\n";
  std::string s1 = "ACCGTT";
  std::string s2 = "ACTTCT";
  size_t ed = editDistance(s1,s2);
  std::cout << ed << "\n\n";
}


// returns the average edit distance between ALL strings with length exactly n
double
testEditDistanceExhaustive(size_t n) {
  double dist = 0;
  char bases[] = {'A','C','G','T'};
  std::string s1(n,'A');
  std::string s2(n,'A');
  s1[1] = 'T';
  for (int i = 0; i < n; ++i) {
    for (int bi = 0; bi < 4; ++bi) {
      s2[i] = bases[bi];
      size_t d = editDistance(s1,s2);
      dist += d;
    }
  }
  return dist;
}

double
recursiveExhEditDistance(std::string s1, std::string s2, size_t n) {
  if (s1.size() == n) {
    double d = editDistance(s1,s2);
    //    std::cout << s1 << '\t' << s2 << '\t' << d << '\n';
    return d;
  } else {
    double ed =
      recursiveExhEditDistance(s1 + "A", s2, n) +
      recursiveExhEditDistance(s1 + "C", s2, n) +
      recursiveExhEditDistance(s1 + "G", s2, n) +
      recursiveExhEditDistance(s1 + "T", s2, n);
    return ed;
  }
}

double recursiveEditDistAllPairs(std::string s1, std::string s2, size_t n) {
  if (s2.size() == n) {
    return recursiveExhEditDistance(s1,s2,n);
  } else {
    double ed =
      recursiveEditDistAllPairs("", s2 + "A", n) +
      recursiveEditDistAllPairs("", s2 + "C", n) +
      recursiveEditDistAllPairs("", s2 + "G", n) +
      recursiveEditDistAllPairs("", s2 + "T", n);
    return ed;
  }
}

void test2bitsEncoding() {
  std::string tests[] = { "ACCTC", "GGG", "", "AAA", "TTTT", "AAGTT", "TTGAA",
			  "AACCCTGTGGACGTGTGACGTGTTGCAAAGCA" };
  for (int  i = 0; i <8; ++i) {
    std::string s = tests[i];
    size_t n = s.size();
    uint64_t e  = string2Encode(s);
    std::string rev_e = encoding2String(e,n);
    std::cout << s << "\t" << e << '\t' << rev_e << '\t' << (rev_e == s) << '\n';
  }
}


void
editDistanceTests() {
  std::cout << "\nEDIT DISTANCE TESTS\n" << std::endl;

  std::cout << "* Reverse and symmetry strings test\n\n";
  for (int i = 0; i < 1000; ++i) {
    size_t n = (rand() % 150 ) + 1;
    size_t m = (rand() % 100 ) + 1;
    std::string s1(n,'N');
    std::string s2(m,'N');
    generateIIDString(s1);
    generateIIDString(s2);

    size_t ed, ed_sym, ed_rev;
    ed = editDistance(s1,s2);
    ed_sym = editDistance(s2,s1);
    std::string s1r(s1.size(), 'N');
    std::string s2r(s2.size(), 'N');
    for (int i = s1.size() - 1, j=0; i >= 0; --i,++j) {
      s1r[j] = s1[i];
    }
    
    for (int i = s2.size() - 1, j=0; i >= 0; --i,++j) {
      s2r[j] = s2[i];
    }
    ed_rev = editDistance(s1r,s2r);

    if ( ed != ed_sym || ed != ed_rev) {
      std::cout << s1 << '\t' << s2 << '\t' << ed << '\t' << ed_sym << '\t' << ed_rev << '\n';
    }
  }

  std::cout << "* Linear space test\n\n";

  size_t v0[200], v1[200];
    
  for (int i = 0; i < 1000; ++i) {
    memset(&v0, 0, 200*sizeof(size_t));
    memset(&v1, 0, 200*sizeof(size_t));
    size_t n = (rand() % 150 ) + 1;
    size_t m = (rand() % 100 ) + 1;
    std::string s1(n,'N');
    std::string s2(m,'N');
    generateIIDString(s1);
    generateIIDString(s2);
    size_t ed, ed_lin;
    ed = editDistance(s1,s2);
    ed_lin = editDistanceLinSpace(s1,s2,v0,v1);
    if (ed != ed_lin) {
      std::cout << s1 << '\t' << s2 << '\t' << ed << '\t' << ed_lin << '\n';
    }
  }
    
  std::cout << "* Edit script test\n\n";
  std::string sa = "ACCTG";
  std::string sb = "GCTGGC";
  EditDistanceInfo info;
  editDistanceWithInfo(sa, sb, info);
  std::cout << info.edit_script << "\t" << info.n_sub << " "
	    << info.n_ins   << " " << info.n_del << "  " 
	    << info.distance() << '\n';
  
  std::cout << "\n* Linear info test\n\n";
  EditDistanceInfo v00[100];
  EditDistanceInfo v11[100];


  /* 	GAGCAACC	TTCGGCGA
  //  editDistanceLinSpaceInfo(std::string("GTCAATGG"), std::string("CCGTTATA"), v00, v11);
  EditDistanceInfo infLin, infQuad;
  infLin = editDistanceLinSpaceInfo(std::string("AGGCCCCT"), std::string("TTCCAATG"), v00, v11);
  std::cout << "\n\n\n";
  editDistanceWithInfo(std::string("AGGCCCCT"), std::string("TTCCAATG"), infQuad);
  std::cout << infLin << std::endl << infQuad << std::endl;*/

  std::string st1(8,'N');
  std::string st2(8,'N');
  for (int i = 0; i < 200; ++i) {
    generateIIDString(st1);
    generateIIDString(st2);
    EditDistanceInfo in = editDistanceLinSpaceInfo(st1, st2, v00, v11);
    EditDistanceInfo infoCorrect;
    editDistanceWithInfo(st1,st2, infoCorrect);
    if (in != infoCorrect) {
      std::cout << st1 << "\t" << st2 << "\n";
      std::cout << in << std::endl << infoCorrect << "\n\n";
    }
  }

  std::cout << "* Bounded Estimate Test" << "\n\n";
  size_t old_opt_k = Options::opts.k;
  Options::opts.k = 3000;
  size_t n = 200;
  double prec = 0.1;
  double z = 2;
  SampleEstimates estimates = editDistanceErrorBoundedEstimates(n, prec, z);
  std::cout << estimates.sampleSize << "\t" << estimates.sampleMean << "\t"
	    << estimates.sampleVariance << std::endl;
  Options::opts.k = old_opt_k;

  std::cout << "\n* Relative Error Estimate Test\n\n";
  old_opt_k = Options::opts.k;
  Options::opts.k = 300000000;
  n = 1000;
  double e_model = 0;

  // n = 35000;
  // double e_model = 18070;
  prec = 0.001;
  z = 2;
  estimates = editDistanceRelativeErrorEstimates(n, e_model, prec, z);
  std::cout << estimates.sampleSize << "\t" << estimates.sampleMean << "\t"
	    << estimates.sampleVariance << "\t" << std::abs(estimates.sampleMean - e_model) << "\n";

  std::cout << "\n* Relative Error Difference Test\n\n";

  prec = 0;
  z=2;
  n = 128;
  std::vector<SampleEstimates> allEst = differenceBoundedRelativeErrorEstimate(n, prec, z, Options::opts.k);
  estimates = allEst[1];

  std::cout << estimates.sampleSize << "\t" << estimates.sampleMean << "\t"
	    << estimates.sampleVariance << "\n";

  Options::opts.k = old_opt_k;


  std::cout << "\n* Sample Matrix Test\n\n";
  
  n = 8;
  EditDistanceInfo** sMat = new EditDistanceInfo*[n];
  for (size_t i = 0; i < n; ++i) {
    sMat[i] = new EditDistanceInfo[n];
  }
  editDistSamplesInfoLinSpace(n,100, sMat);
  //  printMatrix<EditDistanceInfo>(sMat,n,n, "\t");
  for (size_t i = 0; i < n; ++i) {
    delete[] sMat[i];
  }
  delete[] sMat;
}


class InverseSquareRootFunction
{
private:
  double gamma, beta;
public:
  InverseSquareRootFunction(double gam, double bet) {
    gamma = gam;
    beta = bet;
  }
  double operator() (double x) { return ( gamma / std::pow(x, beta) ); }
};

// TODO Move to edit.hpp as soon as good
// This will eventally become the code to verify the model
// E[ed(X,Y)] = \xi n + \gam n^{1-\bet}
void
editDistanceVerifySecondOrderFunction() {
  size_t n_max = std::pow(2,14);
  size_t n = std::pow(2,0);
  InverseSquareRootFunction g(1.0, 2.0);
  double prec = 0.001; //Options::opts.precision;
  double z = 2; //Options::opts.confidence;
  Options::opts.k = 50000;
  while (n <= n_max) {
    SampleEstimates est =
      editDistanceRelativeErrorEstimates(n, 0, prec, z);
    std::cout << n << "\t" << est.sampleMean << "\t" << est.sampleVariance << "\t" << est.sampleSize << std::endl;
    n *= 2;
    
  }
  std::cout << std::endl;
}

void
testAverageDPMatrix(size_t n) {
  double** dpMatrix = new double*[n+1];
  for (int i = 0; i <= n; ++i) {
    dpMatrix[i] = new double[n+1];
  }

  computeAverageDPMatrix(dpMatrix, n, n);

  for (int i = 0; i <= n; ++i) {
    std::cout << dpMatrix[i][i] << std::endl;
  }

  for (int i = 0; i <= n; ++i) {
    for (int j = 0; j <=n; ++j) {
      std::cout << dpMatrix[i][j] << "\t";
    }
    std::cout << "\n";
  }

  for (int i = 0; i <= n; ++i) {
    delete[] dpMatrix[i];
  }
  delete[] dpMatrix;
}

void testAll() {
  std::cout << "--------------------------------\n";
  std::cout << "          TESTING MODE          \n";
  std::cout << "--------------------------------\n";
  //testUtils();
  //testProbFunctions();
  
  //  testSampleEstimators();
  //  testScoreFunction();
  //  testLookupTables();
  //testPeq();
  //testApproximatedExpectedScore();
  
  
  editDistanceTests();

  //editDistanceVerifySecondOrderFunction();
  //  testAverageDPMatrix(Options::opts.N);
}
