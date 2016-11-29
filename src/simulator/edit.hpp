#ifndef SIM_EDIT_H
#define SIM_EDIT_H

#include <memory>
#include <iostream>

// DEPRECATED
#define EDIT_DISTANCE_MONTE_CARLO 1
#define EDIT_DISTANCE_MONTE_CARLO_COMPLETE 2
#define EDIT_DISTANCE_EXHAUSTIVE_ENC 4

#define EDIT_DISTANCE_ESTIMATE_EXHAUSTIVE 0x1
#define EDIT_DISTANCE_ALGORITHM_QUADRATIC 0x2
#define EDIT_DISTANCE_INFO_PARTIAL        0x4
#define EDIT_DISTANCE_INFO_SCRIPT         0x8
#define EDIT_DISTANCE_SAMPLE_MATRIX       0xF


/**
 * This is a structure to store information about how edit distance
 * is divided into substitution, 
 */
class EditDistanceInfo
{
public:
  size_t n_sub;
  size_t n_ins;
  size_t n_del;

  std::string edit_script = "";

  size_t distance() { return n_sub + n_ins + n_del; }
  void reset() { n_sub = 0; n_del = 0; n_ins = 0; }
  bool operator==(const EditDistanceInfo& i) {
    return (n_sub == i.n_sub && n_del == i.n_del && n_ins == i.n_ins);
  }
  bool operator!=(const EditDistanceInfo& i) {
    return !(*this == i);
  }
  
  friend std::ostream& operator<<(std::ostream& out, EditDistanceInfo& info) {
    out << info.n_sub << " " << info.n_del << " " << info.n_ins;
    return out;
  }
  
};

/**
 * \brief computes the edit distance between strings s1 and s2
 */
size_t editDistance(const std::string& s1, const std::string& s2);

size_t editDistanceEncoded(uint64_t s1, size_t n1, uint64_t s2, size_t n2, size_t** dpMatrix);

size_t editDistanceLinSpace(const std::string& s1, const std::string& s2, size_t* v0, size_t* v1);

EditDistanceInfo editDistanceLinSpaceInfo(const std::string& s1, const std::string& s2, EditDistanceInfo* v0, EditDistanceInfo* v1);

/**
 * \brief Computes the dynamic programming matrix dpMatrix for edit distance
 * between s1 and s2
 */
void editDistanceMat(const std::string& s1, const std::string& s2, size_t** dpMatrix);

void editDistanceWithInfo(const std::string& s1, const std::string& s2, EditDistanceInfo& info);

/**
 * \breif Reconstruct the edit script from the dynamic programming matrix.
 * The script and other informations are stored in the passed EditDistanceInfo 
 * structure.
 *
 * \param dpMatrix the dynamic programming matrix
 * \param n number of rows of the matrix minus one 
 * \param m number of columns of the matrix minus one
 * \param info the structure that will be filled with the information
 * 
 */
void
editDistanceBacktrack(size_t** dpMatrix,const std::string& s1, const std::string& s2, EditDistanceInfo& info);


/**
 * Use a Monte-Carlo sampling technique to estimate the edit distance between
 * random strings of same length.
 *
 * \param n_min The minimum length to be tested
 * \param n_max The maximum lengths to be tested
 * \param n_step The incremente on length 
 * \param k_max The number of samples for each length
 * 
 */
void
editDistanceEstimations(size_t n_min, size_t n_max, size_t n_step, size_t k_max);

/*
 * \brief Computes k_the edit distance for k_samples pairs of random strings
 * with length n. The distances are stored in an array whose (smart) pointer 
 * is returned.
 *
 * \param n the length of the strings
 * \param k_samples the number of pairs to compute
 * \return a (smart) pointer to an array containing all the distances
 */
std::unique_ptr<size_t[]>
editDistSamples(size_t n, size_t k_samples);

std::unique_ptr<EditDistanceInfo[]>
editDistSamplesInfo(size_t n, size_t k_samples);


std::unique_ptr<EditDistanceInfo[]>
editDistSamplesInfoLinSpace(size_t n, size_t k_samples, double** avgMatrix = NULL);


double
testExhaustiveEditDistanceEncoded(size_t n);

void
computeAverageDPMatrix(double** dpMatrix, size_t n, size_t m);

#endif
