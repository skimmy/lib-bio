#include "CSFastFormat.hpp"

#define MAX_BUFFER_SIZE 165536

/*

  

*/

/*************** CONSTRUCTORS AND DESTRUCTOR ****************/

CSFastFormat::CSFastFormat(const string& basesFileName, const string& qualFileName) 
  : Format("CSFAST"), sequence(""), header(""), 
    baseFile(basesFileName), qualFile(qualFileName)
{
  nextRead = 0;
  loadNextRead();
}

CSFastFormat::~CSFastFormat() {
  if (nextRead) {
    delete nextRead;
  }
}


/***************** FORMAT METHODS OVERRIDE ******************/

string CSFastFormat::loadFromFile(const string &fileName) {
  ifstream ifs(fileName);
  char buffer[MAX_BUFFER_SIZE];
  sequence.clear();
  while(!ifs.eof()) {
    ifs.getline(buffer, MAX_BUFFER_SIZE);
    sequence += string(buffer);
  }
  return sequence;
}

string CSFastFormat::getSequence() const {
  return sequence;
}

string CSFastFormat::getHeader() const {
  return header;
}

/*********************** LOAD METHODS ***********************/

bool CSFastFormat::hasNextRead() {
  return (nextRead != 0 && nextRead->getBases().length() > 0);
}

CSFastRead CSFastFormat::getNextRead() {
  CSFastRead toReturn = (nextRead) ? CSFastRead(*nextRead) : CSFastRead();
  loadNextRead();
  return toReturn;
}

/********************** UTILITY METHODS *********************/

void CSFastFormat::loadNextRead() {
  if(nextRead) {
    delete nextRead;
  }
  nextRead = new CSFastRead();
  if(baseFile.is_open() && qualFile.is_open()) {
    nextRead->loadBasesAndQualitiesFromFiles(baseFile, qualFile);
  }
}

/************************************************************/
