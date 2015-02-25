/*
 *
 * \file QualityCommon.hpp
 *
 */

#ifndef QUALITY_COMMON_H
#define QUALITY_COMMON_H

#include <string>

namespace qual
{

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

extern const double PHRED[];

const int SANGER_ASCII_OFFSET = 33;
const int SOLEXA_ASCII_OFFSET = 64;
const int ILLUMINA_ASCCI_OFFSET = 64;

enum class QualityEncodingType { SANGER, ILLUMINA, SOLEXA, UNKNOWN };

  /*
   * \fn QualityEncodingType parseQuality(const std::string& quals);
   * \brief Parse the given qualities string and trys to determine which encoding
   * of quality values was used
   */
  QualityEncodingType parseQuality(const std::string& quals);

  /*
   * \fn   std::string encodingToString(QualityEncodingType encoding);
   * \brief Returns a string representation of the quality encoding type
   */
  std::string encodingToString(QualityEncodingType encoding);

}

#endif
