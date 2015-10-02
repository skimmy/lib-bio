#include "../quality.h"

namespace qual
{

  QualityEncodingType parseQuality(const std::string& quals) {
    QualityEncodingType type = QualityEncodingType::SANGER;
    for (size_t i = 0; i < quals.size(); ++i) {
      char q = (char) quals[i]; 
      // all ASCII chars less than 59 are used only by sanger encoding
      if (q < 59) {
	return QualityEncodingType::SANGER;
      }
      // if code is greater than 73 we return Illumina (we don't support Solexa)
      if (q > 73) {
	return QualityEncodingType::ILLUMINA;
      }
      // these values are errors because not covered by any encoding
      if (q < 33 || q > 126) {
	return QualityEncodingType::UNKNOWN;
      }
    }
    return type;
  }

  std::string encodingToString(QualityEncodingType encoding) {
    switch(encoding) {
    case QualityEncodingType::SANGER:
      return "SANGER";
      break;
    case QualityEncodingType::ILLUMINA:
      return "ILLUMINA";
      break;
    case QualityEncodingType::SOLEXA:
      return "SOLEXA";
      break; 
    default:
      return "UNKNOWN";
    }
  }

}
