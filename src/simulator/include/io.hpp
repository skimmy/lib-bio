#ifndef SIM_IO_H
#define SIM_IO_H

#include <string>
#include <list>

typedef std::pair<std::string, int64_t> AlignPair;

std::string loadFromFile(const std::string &fileName);

void loadAlignFromSAM(const std::string &filePath, std::list<AlignPair> &aligns);

#endif
