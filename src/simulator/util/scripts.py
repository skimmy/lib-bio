#!/usr/bin/python

import argparse
import numpy as np

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input scripts file")
    parser.add_argument("n", help="Length of strings")
    parser.add_argument("--matrix-file", dest="matFile", help="Name of file where matrix will be saved", default="matrix.txt")
    parser.add_argument("--max-file", dest="maxFile", help="Name of the file where maximum oscillation vectore will be saved",
                        default="max.txt")
    parser.add_argument("--dist-file", dest="distFile", help="Name od file where distributions of operations will be saved",
                        default="dist.txt")
    # parser.add_argument("opt", help="Required option")
    # parser.add_argument("-o", "--optional", dest='o', help="Optional flag", action='store_true')
    # parser.add_argument("-d", "--default", help="SWith default", default="Hello")
    return parser.parse_args()

def scriptOperationsParse(script):
    matches = 0
    substitutions = 0
    deletions = 0
    insertions = 0
    for op in script:
        if op == 'M':
            matches += 1
        if op == 'S':
            substitutions += 1
        if op == 'D':
            deletions += 1
        if op == 'I':
            insertions += 1
    return (substitutions + deletions + insertions,
            matches, substitutions, deletions, insertions)

def printScripts(scripts):
    for script in scripts:
        (dist,m,s,d,i) = scriptOperationsParse(script)
        print("{0}\t{1}".format(script,dist))

def printFrequencyMatrix(freqMat):
    for row in freqMat:
        print('\t'.join([str(x) for x in row]))

def operationsDistribution(n, scripts):
    subDist = np.zeros(n)
    insDist = np.zeros(n)
    delDist = np.zeros(n)
    matchDist = np.zeros(n)
    for script in scripts:
        nSub = script.count('S')
        nIns = script.count('I')
        nDel = script.count('D')
        nMatch = len(script) - (nSub + nIns + nDel)
        subDist[nSub] += 1
        insDist[nIns] += 1
        delDist[nDel] += 1
        matchDist[nMatch] += 1
    return (matchDist, subDist, delDist, insDist)

def calculateScriptsStats(n, scripts):
    mat = np.zeros((n,n))
    maxDistr = np.zeros(n)
    for s in scripts:
        i = 0
        j = 0
        maxOscill = 0
        mat[i][j] += 1
        for op in s:
            if (op == 'M') or (op == 'S'):
                i += 1
                j += 1
            if (op == 'D'):
                i += 1
            if (op == 'I'):
                j += 1
            mat[i][j] += 1
            if abs(j-i) > maxOscill:
                maxOscill = abs(j-i)
        maxDistr[maxOscill] += 1
    return mat, maxDistr

if __name__ == "__main__":
    args = parseArguments()
    n = int(args.n)
    inFile = args.input
    fh = open(inFile, "r")
    allScripts = []
    for script in fh:
        allScripts.append(script.strip())
    fh.close()
    # n+1 because the DP goes from (0,0) to (n,n) included
    freqMat, maxDistr = calculateScriptsStats(n+1, allScripts)
    #    printScripts(allScripts)
    #    printFrequencyMatrix(freqMat)
    mDist, sDist, dDist, iDist = operationsDistribution(n+1, allScripts)
    np.savetxt(args.distFile ,zip(np.arange(n+1), mDist, sDist, dDist, iDist), delimiter="\t")
    np.savetxt(args.matFile, freqMat, delimiter="\n")
    np.savetxt(args.maxFile, maxDistr, delimiter="\n")
    
    
    
    
