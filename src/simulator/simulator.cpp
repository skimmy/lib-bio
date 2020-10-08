/**
 * This is a simulator developed ad-hoc to allow future optimization
 * regardless the changes to other components of the library (which
 * are designed ti be part of a library rather then efficient stand
 * alon tools.
 */

#include <include/common.hpp>
#include <include/options.hpp>

// Basic includes
#include <include/generator.hpp>
#include <include/log.hpp>
#include <include/prob.hpp>
#include <include/util.hpp>


// Task includes
#include <include/align.hpp>
#include <include/chain.hpp>
#include <include/edit.hpp>
#include <include/edit_estimates.hpp>
#include <include/online.hpp>
#include <include/lcs.hpp>

// standard includes
#include <cstdlib>
#include <ctime>

#include <unordered_map>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>

// TODO Convert to use proper ns
using namespace lbio::sim;
using namespace lbio::sim::generator;
using namespace lbio::sim::edit;
using namespace lbio::sim::lcs;
using namespace lbio::sim::log;

// output quantities (common to online and offline)
double p_fail = 0.0;
size_t holes = 0;
size_t actually_produced_reads = 0;
double scoreSum = 0.0;
EmpiricalDistribution scoreDist(0, 1, 10);

// output quantities for oracle simulations

/* 
   An estimation point contains for each s=0,...,m
     - sum of all scores
     - sum of all denominators
     - sum of all denominators
     - number of observations of 's'
   This can then be used to compute average quantities and compare
   them (i.e., esitmated) with the theoretical one (possibly
   approximated).
*/
struct EstimationPoint {
    double sumNum;
    double sumDen;
    double sumScore;
    int count;
    int hammDist;
};
EstimationPoint *oraclePoints;

EditDistanceSimOutput *edOut;

void initSimulator() {
    initUtil(Options::opts.m);
    initRandomGenerator(Options::opts.alphabet.size(), Options::opts.seed);
    initProbabilities();
    initChainMatrix();

    scoreDist
            = EmpiricalDistribution(0, 1, Options::opts.empiricalDistributionStep);
    oraclePoints = new EstimationPoint[Options::opts.m + 1];

    for (size_t i = 0; i < Options::opts.m; ++i) {
        oraclePoints[i].sumNum = 0.0;
        oraclePoints[i].sumDen = 0.0;
        oraclePoints[i].sumScore = 0.0;
        oraclePoints[i].count = 0;
        oraclePoints[i].hammDist = 0;
    }

    // set the precision from options
    std::cout.precision(Options::opts.floatPrecision);

    edOut = new EditDistanceSimOutput();
}

void clearSimulator() {
    delete edOut;
    delete[] oraclePoints;
    clearChainMatrix();
    clearProbabilities();
    clearUtil();
}

void outputResults() {
    if (Options::opts.task == Task::EditDist) {
        // CDF ouput requested (-D <distfile> option)
        if (!Options::opts.outputDistribution.empty() &&
            Options::opts.subTask != EDIT_DISTANCE_SUBTASK_SCRIPT_DIST) {
            if (edOut->distPDF) {
                size_t n = Options::opts.N;
                std::ofstream ofs(Options::opts.outputDistribution, std::ofstream::out);
                logInfo("Writing Edit Distance distribution on "
                        + Options::opts.outputDistribution);
                for (size_t i = 0; i <= n; ++i) {
                    ofs << edOut->distPDF[i] << std::endl;
                }
                ofs.close();
            }
        }
        return;
    }
    if (Options::opts.pipeline) {
        if (Options::opts.task == Task::Oracle) {
            double appNumDen[2];
            for (size_t i = 0; i < Options::opts.m + 1; ++i) {
                approximatedScore(i, appNumDen);
                std::cout << i << "\t" << oraclePoints[i].sumScore << "\t"
                          << oraclePoints[i].count << "\t" << oraclePoints[i].sumNum
                          << "\t" << oraclePoints[i].sumDen << "\t" << appNumDen[0]
                          << "\t" << appNumDen[1] << "\t" << oraclePoints[i].hammDist
                          << "\n";
            }
            return;
        }
        std::cout << p_fail << std::endl;
    } else {
        std::cout << "P[Fail]    = " << p_fail << std::endl;
        std::cout << "P[Success] = " << 1.0 - p_fail << std::endl;
        std::cout << "#[Holes]   = " << holes << std::endl;
        std::cout << "#[Reads]   = " << actually_produced_reads << std::endl;
    }
    if (!Options::opts.outputDistribution.empty()) {
        std::ofstream ofs(Options::opts.outputDistribution, std::ofstream::out);
        for (size_t i = 0; i < scoreDist.getIntervalCount(); ++i) {
            ofs << scoreDist.valueAtIndex(i) << '\n';
        }
        ofs.close();
    }

    if (!Options::opts.outputCDF.empty()) {
        std::ofstream ofs(Options::opts.outputCDF, std::ofstream::out);
        std::vector<double> cdf(scoreDist.getIntervalCount());
        scoreDist.getCDF(cdf);
        for (size_t i = 0; i < scoreDist.getIntervalCount(); ++i) {
            ofs << cdf[i] << "\n";
        }
        ofs.close();
        std::cout << scoreDist.valueAtIndex(percentileIndex(cdf, 0.001)) << "\n";
    }
}

void recordScoreWithOverlap(double sc, size_t s) {
}

void recordScore(double p_ab) {
    p_fail += 1.0 - p_ab;
    scoreSum += p_ab;
    scoreDist.addSample(p_ab);
}

void offlineSimulation() {

    char *ref = NULL;

    size_t N = Options::opts.N;
    size_t m = Options::opts.m;

    ref = new char[N];
    generateIIDGenome(N, ref);
    std::string s(ref);

    // priority queue is used with position as key so that while
    // extractin reads at once (i.e., emptying the queue) reads will be
    // presented in ordered by position on the reference sequence
    std::priority_queue<Read> reads;
    generateOfflineReads(s, reads, m, Options::opts.M, Options::opts.pe);

    // Temporary variables to count the number of holes, in the future a
    // more sophisticated way (e.g., finite state machine) should be
    // used.
    bool onHole = false;

    Read r1 = reads.top();
    reads.pop();
    while (!reads.empty()) {
        Read r2 = reads.top();
        reads.pop();
        size_t s = m - (r2.j - r1.j);
        evaluateChainRelation(r1, r2, s);

        if (s <= m) {
            onHole = false;
            double p_ab = randomReadsOverlapProbNoErr(r1.r, r2.r, s);
            recordScore(p_ab);

        } else {
            if (onHole == false) {
                onHole = true;
                holes++;
            }
            addNonOverlapRecord(r2.j - r1.j - m);
            double x = (double) Options::opts.N - 2.0 * (double) Options::opts.m + 1.0
                       + overlappingStringsSum(r1.r, r2.r);
            recordScore(1.0 / x);
        }
        r1 = r2;
    }

    delete[] ref;
}

void onlineSimulation() {

    bool onHole = false;

    size_t N = Options::opts.N;
    size_t m = Options::opts.m;
    double pe = Options::opts.pe;

    GenomeSegment g(N, m, MAX_GENOME_SEGMENT_LENGTH);
    generateFirstGenomeSegment(g);

    size_t generated_reads = 0;
    size_t current_position = 0;
    size_t real_position = 0;
    size_t remaining_genome = g.length;

    size_t actual_M = 0;

    Read prev_read("", -1);

    while (real_position < N - m) {

        size_t d = generateInterReadDistance();
        real_position += d;

        // this is artificial however for reasonable values of parameters
        // it should never happen otherwise we woul need a different way
        // of online generating the genome.  More specifically if that
        // happens it means that 'd' is higher then a whole genome segment
        // (which should be no less than 10000 in practical cases) for
        // reasonable values of N and M this event will have probability
        // zero for all practical situations and artifically skipping over
        // such 'extreme' values of d will not appreciably change final
        // results
        if (d > (g.length - m - 1)) {
            continue;
        }

        // in this case we need to generate new genome segment...
        if (remaining_genome < m + d) {
            size_t tmp = current_position + d;
            if (tmp < g.length) {
                generateNewGenomeSegment(g, g.length - tmp);
            } else {
                generateNewGenomeSegment(g, 0);
            }
            current_position = 0;
            remaining_genome = g.length;
        } else {
            // ...otherwise we simply update counters
            current_position += d;
            remaining_genome -= d;
        }

        Read current = generateOnlineRead(g.genome, current_position, m, pe);
        actual_M++;
        current.j = real_position;

        // here the probabilities are computed and accumulated
        if (prev_read.j != (size_t) -1) {

            if (d > m) {
                if (!onHole) {
                    holes++;
                }
                onHole = true;
                // non-overlap case...
                if (Options::opts.approxLevel < 0) {
                    double sc = score(prev_read.r, current.r, 0);
                    recordScore(sc);
                }

            } else {
                onHole = false;
                // overlap case...
                size_t s = m - d;
                double sc = score(prev_read.r, current.r, s);
                recordScore(sc);
            }
        }

        generated_reads++;
        current.j = real_position;
        prev_read = current;
    }
    actually_produced_reads = actual_M;
}

void oracleSimulation() {
    size_t n = 2 * Options::opts.m;
    size_t m = Options::opts.m;
    double pe = Options::opts.pe;
    double alpha = 1.0 / ((double) Options::opts.N - 2.0 * m + 1);
    double numDen[2];

    char *genome = new char[n];
    // Oracle simulation loops to produce exactly M-1 consecutive pairs
    for (size_t i = 0; i < Options::opts.M - 1; ++i) {
        // generate 2m bases of genome
        generateIIDGenome(n, genome);

        // generate first reads at position 0
        Read r1 = generateOnlineRead(genome, 0, m, pe);
        // generate inter-arrival d
        size_t d = generateInterReadDistance();

        // generate second reads at position d
        if (d >= m) {
            oraclePoints[0].sumNum += 1.0;
            oraclePoints[0].sumDen += (1.0 / alpha);
            oraclePoints[0].sumScore += alpha;
            oraclePoints[0].count++;
        } else {
            Read r2 = generateOnlineRead(genome, d, m, pe);
            size_t s = m - d;
            oraclePoints[s].hammDist += prefixSuffixHammingDistance(r1.r, r2.r, s);
            double sc = scoreExt(r1.r, r2.r, s, numDen);
            oraclePoints[s].sumNum = numDen[1];
            oraclePoints[s].sumDen = numDen[0];
            oraclePoints[s].sumScore += sc;
            oraclePoints[s].count++;

        }
    }
    delete[] genome;
}

void
editDistanceOpMode() {


    using AlgorithmBand = EditDistanceBandApproxLinSpace<lbio_size_t, std::string>;
    using AlgorithmExact = EditDistanceWF<lbio_size_t, std::string>;


    // The default edit distance mode is
    // Sample
    // Linear Alg
    // mean and variance output
    // no script
    // no sample matrix
    int flags = Options::opts.optFlags; // -f
    size_t n = Options::opts.N;
    int task = Options::opts.subTask; // -B
    edOut->distPDF = new double[n + 1];
    std::fill_n(edOut->distPDF, n + 1, 0);

    std::vector<lbio_size_t> weights = {1, 1, 1};
    if (Options::opts.weight_scheme == 1) {
        weights = {2 * (n + 1), 1, 1};
    }

    // TASK - Scripts Generation (8)
    if (task == EDIT_DISTANCE_SUBTASK_SCRIPT_DIST) { // -B 8
        logInfo("Task 'Script Distribution'");

        AlgorithmExact alg{n, n, weights};
        std::vector<std::string> allScripts{};
        std::vector<std::pair<std::string, std::string>> allStrings{};

        // verbose > 0 --> more info
        if (Options::opts.verbose > 0) {
            generate_scripts_with_strings(n, n, Options::opts.k, allScripts,
                                          allStrings, alg, Options::opts.alphabet);
            // TODO save to file

            // Or on the standard out
            std::cout << "x,y,s,d,d1,d2\n";
            for (lbio_size_t i = 0; i < Options::opts.k; ++i) {
                std::string x = allStrings[i].first;
                std::string y = allStrings[i].second;
                std::string s = allScripts[i];
                // first and second half distances
                lbio_size_t d = alg.calculate(x, y);
                lbio_size_t d1 = alg.calculate(x.substr(0, n / 2), y.substr(0, n / 2));
                lbio_size_t d2 = alg.calculate(x.substr(n / 2, n / 2), y.substr(n / 2, n / 2));
                std::cout << x << "," << y << "," << s << "," << d << "," << d1 << "," << d2 << "\n";
            }
        } else { // verbose = 0 --> only scripts
            generate_scripts(n, n, Options::opts.k, allScripts, alg, Options::opts.alphabet);
            // if file is given save there otherwise use std out
            if (!Options::opts.outputDistribution.empty()) {
                std::ofstream ofs(Options::opts.outputDistribution, std::ofstream::out);
                for (std::string script : allScripts) {
                    ofs << script << "\n";
                }
                ofs.close();
            } else {
                for (std::string script : allScripts) {
                    std::cout << script << "\n";
                }
            }
        }
        return;
    }

    // TASK - Algorithms comparison (32)
    if (task == EDIT_DISTANCE_SUBTASK_COMPARE_ALGS) {
        logInfo("Task 'Algorithms comparison'");
        compare_edit_distance_algorithms(n, n, Options::opts.k, Options::opts.alphabet);
        return;
    }

    // TASK - Edit distance calculation (0, default)
    if (flags & EDIT_DISTANCE_DIFF_BOUNDED_ERROR) { // -f 64
        logInfo("Task 'g(n) Esitmation'");
        size_t k_max = Options::opts.k;
        double precision = Options::opts.precision;
        double z_confidence = Options::opts.confidence;
        // lambda for output (opt_str for extra algorithm output)
        auto print_cb =
                [=](const SampleEstimates &est_n, const SampleEstimates &est_n_2,
                    std::string opt_str = "") {
                    std::cout << (n >> 1) << "\t" << est_n_2 << "\n"
                              << n << "\t" << est_n << "\t" << opt_str << "\n";
                };

        lbio_size_t T = static_cast<lbio_size_t>(std::floor(n / 2.0));
        lbio_size_t Tmin = static_cast<lbio_size_t>(std::sqrt(n));
        // Approximation is required find 'optimal' T >= sqrt(n)
        if (Options::opts.approxLevel == 1) { // -A 1 --> Band with optimal width estimation
            if (flags & EDIT_DISTANCE_BANDWIDTH_ESTIMATE) {
                logInfo("Estimation of optimal bandwidth...");
                Tmin = std::max(1, Options::opts.approxLevel);
                T = optimal_bandwidth(n, Options::opts.alphabet, precision / 2, 16, Tmin);
                logDebug(debug_string("[Band Estimate]\t", "~T*: " + std::to_string(T)));
            } else {
                T = Options::opts.approxLevel;
            }
            AlgorithmBand alg(n, n, T, weights);
            logInfo("Estimation...");
            std::vector<SampleEstimates> est =
                    edit::difference_estimate(n, precision, z_confidence,
                                              k_max, alg, print_cb, Options::opts.alphabet);
        } // -A 1

        if (Options::opts.approxLevel == 2) { // -A 2 Adaptive T
            logInfo("Adaptive Estimation...");
            std::vector<SampleEstimates> est =
                    edit::difference_estimate_adaptive(n, precision, z_confidence,
                                                       k_max, print_cb, Options::opts.alphabet);
        }
        return;
    }

    if (flags & EDIT_DISTANCE_BOUNDED_ERROR) { // -f 32
        double precision = Options::opts.precision;
        double z_confidence = Options::opts.confidence;
        SampleEstimates beEst =
                editDistanceErrorBoundedEstimates(n, Options::opts.alphabet, precision, z_confidence);
        std::cout << beEst.sampleSize << "\t" << beEst.sampleMean << "\t"
                  << beEst.sampleVariance << "\n";
        return;
    }


    // subtask eccentricity
    if (task == EDIT_DISTANCE_SUBTASK_ECCENTRICITY) { // -B 4

        if (flags & EDIT_DISTANCE_FLAG_ECCENTRICITY_ALL) { // -f 1
            logInfo("Eccentricity all strings");
            lbio_size_t n = Options::opts.N;
            // iterator for all the strings
            AlphabetIterator it(n, Options::opts.alphabet);
            edit_distance_eccentricity(it, it.end(), std::cout, Options::opts.alphabet);
            return;
        }
        if (flags & EDIT_DISTANCE_FLAG_ECCENTRICITY_FILE) { // -f 3
            logInfo("Eccentricity strings from file");
            std::vector<std::string> v;
            std::ifstream ifs(Options::opts.inputReference);
            insert_from_stream(std::inserter(v, v.begin()), ifs);
            edit_distance_eccentricity(v.begin(), v.end(), std::cout, Options::opts.alphabet);
            for (auto it = v.begin(); it != v.end(); ++it) {
                std::cout << eccentricity_for_string(*it, Options::opts.alphabet) << "\n";
            }
            return;
        }
    }


    // Task Exhaustive
    if (flags & EDIT_DISTANCE_ESTIMATE_EXHAUSTIVE) {
        lbio_size_t n = Options::opts.N;
        std::string alphabet = Options::opts.alphabet;
        lbio_size_t t = Options::opts.n_threads;
        EccentricityResult res = eccentricity_with_symmetries(n, alphabet, t);
        double e = std::get<0>(res);
        std::cout << n << "," << e << "," << e / n;
        if (Options::opts.verbose) {
            std::cout << "," << std::get<1>(res) << "," << std::get<2>(res)
                      << "," << std::get<3>(res) << "," << std::get<4>(res);
        }
        std::cout << "\n";
    } else { // Not exhaustive
        // SAMPLE
        size_t k = Options::opts.k;

        if (flags & EDIT_DISTANCE_ALGORITHM_QUADRATIC) {
            // QUADRATIC + Sample
            if (flags & EDIT_DISTANCE_INFO_PARTIAL) {
                // Sample + Quadratic + Partial Info
                if (flags & EDIT_DISTANCE_INFO_SCRIPT) { // -f 14
                    logInfo("Sample quadratic algorithms info");
                    AlgorithmExact algExact{n, n, {1, 1, 1}};
                    EditDistanceSample<AlgorithmExact> generator{n, n, Options::opts.alphabet};
                    for (size_t i = 0; i < k; ++i) {
                        generator(algExact);
                        auto info = algExact.backtrack();
                        std::cout << info.n_sub << "\t" << info.n_ins
                                  << "\t" << info.n_del << "\n";
                    }
                } else {
                    logWarning("Sample quadratic info without script not available");
                }
            }
        } else {
            // LINEAR + Sample

            if (flags & EDIT_DISTANCE_INFO_PARTIAL) { // -f 4
                // PARTIAL INFO + Sample + Linear
                logWarning("Partial info for linear under developement");
                std::unique_ptr<EditDistanceInfo[]> samples =
                        editDistSamplesInfoLinSpace(n, k, Options::opts.alphabet);

                // auto -> std::unique_ptr<double[]>
                auto subSamples = extractSubstitutionArray(samples.get(), k);
                auto delSamples = extractDeletionArray(samples.get(), k);
                auto insSamples = extractInsertionArray(samples.get(), k);

                SampleEstimates subEst
                        = estimatesFromSamples<double>(subSamples.get(), k);
                SampleEstimates delEst
                        = estimatesFromSamples<double>(delSamples.get(), k);
                SampleEstimates insEst
                        = estimatesFromSamples<double>(insSamples.get(), k);

                // If 'verbose' is set  all samples are printed
                if (Options::opts.verbose) {
                    for (size_t i = 0; i < k; ++i) {
                        std::cout << samples[i] << std::endl;
                    }
                }
                std::cout << subEst.sampleMean << "\t" << subEst.sampleVariance << "\n";
                std::cout << delEst.sampleMean << "\t" << delEst.sampleVariance << "\n";
                std::cout << insEst.sampleMean << "\t" << insEst.sampleVariance << "\n";
            }

            if (flags & EDIT_DISTANCE_MINMAX_STRING) { // -f 16
                logInfo("Sampling periodic/constant string distribution");

                std::cout << std::endl;
                // Sampling e("ACGT^*, Y)
                std::string s(n, 'N');
                for (lbio_size_t i = 0; i < n; ++i) {
                    s[i] = bases[i % 4];
                }
                std::vector<size_t> v = edit_samples_fixed_string(n, k, s, Options::opts.alphabet);
                auto est = estimatesFromSamples(v.cbegin(), v.cend(), k);
                std::cout << "ACGT*\t" << est << "\n";

                // sampling  e("A*", Y)
                s = std::string(n, 'A');
                v = edit_samples_fixed_string(n, k, s, Options::opts.alphabet);
                est = estimatesFromSamples(v.cbegin(), v.cend(), k);
                std::cout << "A*\t" << est << "\n";
                std::cout << std::endl;
                return;
            } else { // -f 0
                // MINIMAL INFO (mean + var) + Sample + Linear
                logInfo("Basic sampling");
                lbio_size_t bandwidth = std::ceil(n / 2.0);
                if (Options::opts.approxLevel == 1) { // -A 1 bandwise with sqrt(n) bandwidth
                    bandwidth = std::ceil(std::sqrt(n) / 2.0);
                }
                if (Options::opts.approxLevel == 2) { // -A 2 bandwise with log(n) bandwidth
                    bandwidth = std::ceil(std::log2(n) / 2.0);
                }
                AlgorithmBand alg(n, n, bandwidth, weights);

                std::vector<size_t> v;
                auto ins_ = std::inserter(v, v.begin());
                edit::edit_distance_samples(n, k, ins_, alg, Options::opts.alphabet);
                // save all the info
                if (Options::opts.verbose > 0) {
                    // if file is not indicated print on cerr
                    std::string file_name = Options::opts.verboseOutput;
                    if (file_name == "") {
                        std::copy(v.begin(), v.end(), std::ostream_iterator<size_t>(std::cerr, "\n"));
                    } else {
                        std::ofstream os(file_name);
                        std::copy(v.begin(), v.end(), std::ostream_iterator<size_t>(os, "\n"));
                    }
                }
                auto estimators = estimatesFromSamples(v.begin(), v.end(), k);
                std::cout << estimators.sampleMean << ",";
                std::cout << estimators.sampleVariance << ",";
                std::cout << estimators.sampleThirdAbsMoment << std::endl;
            }
        }
    }
}

void
prototyping() {
    // !!! WARNING: possibly don't remove next two lines !!!
    std::string proto_task_msg = make_bold("");
    logInfo("Working on " + proto_task_msg);
    std::cout << sizeof(lbio_size_t) << "\n";
}

int main(int argc, char **argv) {
    // Important NOT invert (init requires argument to be parsed)


#ifdef HAVE_BOOST_PROGRAM_OPTIONS
    parseArgumentsBoost(argc, argv);
#else
    parseArguments(argc,argv);
#endif


    initSimulator();

    switch (Options::opts.task) {
        case (Task::Test):
            logWarning("Prototyping, for tests run proper binary");
            prototyping();
            exit(0);
        case (Task::Offline):
            offlineSimulation();
            break;
        case (Task::Online):
            onlineSimulation();
            break;
        case (Task::Oracle):
            oracleSimulation();
            break;
        case (Task::AlignScore):
            evaluateAlignmentScore(Options::opts);
            break;
        case (Task::EditDist): // -O 6
            editDistanceOpMode();
            break;
        default:
            std::cout << "Unrecognized operation mode " <<
                      static_cast<int>(Options::opts.task) << "\nAborting..\n";
            exit(1);
    }
    outputResults();
    clearSimulator();

    return 0;
}
