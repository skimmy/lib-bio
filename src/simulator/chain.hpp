#ifndef CHAIN_H
#define CHAIN_H

void computeExpectedFalsePositiveMatrix(double** m);

void printFalsePositiveMatrix();

void initChainMatrix();
void clearChainMatrix();
void printChainMatrix();
void printNonOverlapDistribution();
void evaluateChainRelation(const Read& r1, const Read& r2, size_t s);
void addNonOverlapRecord(size_t d);


#endif
