#ifndef SIM_EDIT_H
#define SIM_EDIT_H

#include <memory>

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

#endif
