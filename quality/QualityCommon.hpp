/*
 *
 * \file QualityCommon.hpp
 *
 */

#ifndef QUALITY_COMMON_H
#define QUALITY_COMMON_H

#include <string>

#define SANGER_ASCII_OFFSET 33
#define SOLEXA_ASCII_OFFSET 64
#define ILLUMINA_ASCCI_OFFSET 64

namespace qual
{

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
