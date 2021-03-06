#include "common.hpp"

#include "generator.hpp"
#include "options.hpp"
#include "prob.hpp"
#include "util.hpp"
#include "edit.hpp"
#include "log.hpp"
#include "edit_estimates.hpp"

#include "extensions/boost_ext/boost_edit.hpp"

#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>

using namespace lbio::sim;
using namespace lbio::sim::generator;

#define TEST_GENOME_LENGTH 65536
#define TEST_READ_LENGTH 100

const int MC_SAMPLES = 2 << 22;

void test_edit_distance_class();

void
testUtils() {
    std::cout << "\nUTILITY FUNCTIONS/CLASSES TESTS" << std::endl;

    std::cout << "* Geometric progression\n";
    GeometricProgression<int> geom(2, 1);
    for (size_t i = 0; i < 10; ++i) {
        std::cout << geom.getNext() << std::endl;
    }

    std::cout << "* Linear Progression\n";
    LinearProgression<double> lin(0.5, 0);
    for (size_t i = 0; i < 10; ++i) {
        std::cout << lin.getNext() << std::endl;
    }
}

void
testProbFunctions() {
    std::cout << "\nPROBABILITY FUNCTIONS TESTS\n\n";

    std::cout << "* Median tests\n";
    int f_sym[5] = {1, 3, 4, 3, 1};
    size_t median_idx = medianFromFrequency<int>(f_sym, 5);
    std::cout << "  Median for {1, 3, 4, 3, 1}\t\t"
              << f_sym[median_idx] << " i: " << median_idx << "\n";

    double f_skew[7] = {0.5, 0.5, 0.5, 1, 1, 2, 0.5};
    median_idx = medianFromFrequency<double>(f_skew, 7);
    std::cout << "  Median for {0.5, 0.5, 0.5, 1, 1, 2, 0.5}\t\t"
              << f_skew[median_idx] << " i: " << median_idx << "\n";

    double f_0[3] = {0.5, 0.25, 0.25};
    median_idx = medianFromFrequency<double>(f_0, 3);
    std::cout << "  Median for {0.5, 0.25, 0.25}\t\t"
              << f_0[median_idx] << " i: " << median_idx << "\n";

}

void
testSampleEstimators() {

    std::cout << "\nESTIMATORS TESTS\n\n";

    std::cout << "Simple test\n Samples: [1, 2, -2, 4, 2, 2]\n";
    size_t k = 6;
    int x[] = {1, 2, -2, 4, 2, 2};
    SampleEstimates est = estimatesFromSamples<int>(x, k);
    std::cout << "* Sample Mean:   " << est.sampleMean << std::endl;
    std::cout << "* Sample Var:    " << est.sampleVariance << std::endl;

    std::cout << "\nConstant test\n Samples = [0.5, 0.5, 0.5]\n";
    double y[] = {0.5, 0.5, 0.5};
    k = 3;
    est = estimatesFromSamples<double>(y, k);
    std::cout << "* Sample Mean:   " << est.sampleMean << std::endl;
    std::cout << "* Sample Var:    " << est.sampleVariance << std::endl;

    int z[128];
    for (int i = 0; i < 128; ++i) {
        z[i] = (rand() % 65) - 32;
    }
    std::cout << "\nRandom test\n 128 samples U[-32,32]\n";
    est = estimatesFromSamples<int>(z, 128);
    std::cout << "* Sample Mean:   " << est.sampleMean << std::endl;
    std::cout << "* Sample Var:    " << est.sampleVariance << std::endl;

}

void testApproximatedExpectedScore() {
    for (int s = 0; s < (int) Options::opts.m + 1; ++s) {
        std::cout << s << "\t" << approximatedScore(s) << "\n";
    }
}

// tests the value of p_eq by simulating 
void testPeq() {

    double pe = Options::opts.pe;
    double jointProb[4][4];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            jointProb[i][j] = 0;
        }
    }
    // joint probability distribution (R1 = x, R2 = y | R0 = tBase)
    double x = 0;
    char R1, R2;
    for (int i = 0; i < MC_SAMPLES; ++i) {
        char tBase = bases[rand() & 0x03];
        R1 = R2 = tBase;
        x = (double) rand() / RAND_MAX;
        if (x < pe) {
            R1 = randomMutation(tBase);
        }
        x = (double) rand() / RAND_MAX;
        if (x < pe) {
            R2 = randomMutation(tBase);
        }
        jointProb[(int) revBases[(int) R1]][(int) revBases[(int) R2]]++;
    }
    std::cout << "\n";
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            std::cout << jointProb[i][j] / (double) MC_SAMPLES << '\t';
        }
        std::cout << '\n';
    }
}


void testLookupTables() {
    clearUtil();
    initUtil(Options::opts.m);
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

    simulateReadAt(0, TEST_READ_LENGTH, S, r1);
    simulateReadAt(0, TEST_READ_LENGTH, S, r2);

    for (size_t i = 1; i <= TEST_READ_LENGTH; ++i) {
        std::cout << i << '\t' << score(std::string(r1), std::string(r2), i) << std::endl;
    }


    // TODO this test is 'trivial' leave as last one
    // Test (ii) two reads with no prefix suffix overlap
    std::cout << "(ii)\n";

    Options::opts.pe = 0;
    Options::opts.m = TEST_READ_LENGTH;
    Options::opts.N = TEST_GENOME_LENGTH;

    simulateReadAt(0, TEST_READ_LENGTH, S, r1);
    simulateReadAt(0, TEST_READ_LENGTH, S, r2);


    // Test (iii) two reads with exactly one prefix suffix overlap
    std::cout << "(iii)\n";

    size_t m = TEST_READ_LENGTH;
    simulateReadAt(0, TEST_READ_LENGTH, S, r1);
    r1[m - 1] = baseComplement(r1[m - 1]);
    r1[0] = 'G';
    // test all overlaps except exact one
    for (size_t s = 1; s < m; ++s) {
        simulateReadAt(0, TEST_READ_LENGTH, S, r2);
        r2[s - 1] = r1[m - 1];
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

void
editDistanceTests() {

}

void testAll() {

    logInfo("TEST MODE");
    //testUtils();
    //testProbFunctions();

    //  testSampleEstimators();
    //  testScoreFunction();
    //  testLookupTables();
    //testPeq();
    //testApproximatedExpectedScore();


    editDistanceTests();
    test_edit_distance_class();
}
