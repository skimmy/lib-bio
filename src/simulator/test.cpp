#include "common.hpp"

#include <iostream>

#define TEST_GENOME_LENGTH 65536
#define TEST_READ_LENGTH 100

const int MC_SAMPLES = 2 << 22;

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

void
editDistanceEstimations(size_t n_min, size_t n_max, size_t n_step, size_t k_max) {
  std::cout << std::endl;
  for (size_t n = n_min; n <= n_max; n += n_step) {
    std::string s1(n,'N');
    std::string s2(n,'N');
    double AED = 0;
    for (size_t k = 1; k <= k_max; ++k) {
      generateIIDString(s1);
      generateIIDString(s2);
      AED += editDistance(s1,s2);
    }
    std::cout << n << "\t" << ( AED / k_max) << std::endl;
  }
  std::cout << std::endl;
}

void testAll() {
  std::cout << "--------------------------------\n";
  std::cout << "          TESTING MODE          \n";
  std::cout << "--------------------------------\n";
  //  testScoreFunction();
  //  testLookupTables();
  //testPeq();
  //testApproximatedExpectedScore();
  editDistanceEstimations(250,10000,25,1000);
}
