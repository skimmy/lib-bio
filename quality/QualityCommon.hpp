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
