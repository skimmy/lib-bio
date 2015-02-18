#ifndef QUALITY_CONST
#define QUALITY_CONST

/**
 *
 * \file phred_const.hpp
 * \brief Contains constant related to <i>phred function</i> and related quality
 * values.
 *
 */

/**
 * \var double PHRED;
 * \brief Pre-computed constants for <i>phred</i> qualities.
 *
 * This array contains pre-computed values of the PHRED function defined by
 * <CENTER>
 *   \f$ PHRED[q] = 10^{-q / 10} \f$
 * </CENTER>
 * This arrya is computed for quality values in the closed range [0,300)
 */

extern double PHRED[];



#endif
