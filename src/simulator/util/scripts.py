#!/usr/bin/python

import argparse
import numpy as np
import os.path as path
import os
import sys

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

redErrorMsg = bcolors.FAIL + "Error" + bcolors.ENDC
yellowWarnMsg = bcolors.WARNING + "Warning" + bcolors.ENDC

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input scripts file")
    parser.add_argument("n", help="Length of strings")
    parser.add_argument("--matrix-file", dest="matFile", help="Name of file where matrix will be saved", default="matrix.txt")
    parser.add_argument("--max-file", dest="maxFile", help="Name of the file where maximum oscillation vectore will be saved",
                        default="max.txt")
    parser.add_argument("--dist-file", dest="distFile", help="Name od file where distributions of operations will be saved",
                        default="dist.txt")
    parser.add_argument("--mat-row", dest="matrixRow", help="Prints the content of the specified row of the frequency matrix")
    parser.add_argument("--mat-col", dest="matrixColumn", help="Prints the content of the specified column of the frequency matrix")
    parser.add_argument("--stretch-stat", dest="stretchStat", default=None,
                        help="Name of the files (one per operation) where statistic of stretches will be saved")
    parser.add_argument("--diagonal-stat", dest="diagonalStats", default=None,
                        help="Name of the files the will contain passages through diagonal statistics of scripts")
    parser.add_argument("-a", "--archive-dir", dest="archive", help="Creates an 'archived' version in the directory given as parameter")
    
    
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

def computeStretchesDistribution(n, scripts):
    dictOfDists = {}
    mSDist = np.zeros(n)
    sSDist = np.zeros(n)
    dSDist = np.zeros(n)
    iSDist = np.zeros(n)
    dictOfDists['M'] = mSDist
    dictOfDists['S'] = sSDist
    dictOfDists['D'] = dSDist
    dictOfDists['I'] = iSDist
    for script in scripts:
        currentChar = script[0]
        countedChar = 0
        for s in script:            
            if (s == currentChar):
                countedChar += 1
            else:
                dictOfDists[currentChar][countedChar] += 1
                currentChar = s
                countedChar = 1
        dictOfDists[currentChar][countedChar] += 1
    return (mSDist, sSDist, dSDist, iSDist)

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

def horizontalCross(script, r):
    i = 0
    j = 0
    for op in script:
        next_i = i
        next_j = j
        if (op == 'M') or (op == 'S') or (op == 'D'):
            next_i = i+1
        if (op == 'M') or (op == 'S') or (op == 'I'):
            next_j = j+1
        if (next_j >= r):
            return (i,j)
        i = next_i
        j = next_j

def calculateHorizontalFlow(n, scripts, r):
    f = np.zeros(n)
    for s in scripts:
        i,j = horizontalCross(s,r)
        f[i] += 1
    return f

def verticalCross(script, c):
    i = 0
    j = 0
    for op in script:
        next_i = i
        next_j = j
        if (op == 'M') or (op == 'S') or (op == 'D'):
            next_i = i+1
        if (op == 'M') or (op == 'S') or (op == 'I'):
            next_j = j+1
        if (next_i >= c):
            return (i,j)
        i = next_i
        j = next_j

def calculateVerticalFlow(n, scripts, c):
    f = np.zeros(n)
    for s in scripts:
        i,j = verticalCross(s, c)
        f[j] += 1
    return f

'''Computes the statistics of scripts with respect to the main diagonal.
More precisely computes the statistics of entering and exiting point (i.e.,
which is the frequency of (i,i) being enter/exit point and the frequency of
lengths of subpaths staying on the diagonal.
'''
def calculateDiagonalStats(n, scripts, diag=0):
    enterStats = np.zeros(n) 
    exitStats = np.zeros(n)
    onDiagStats = np.zeros(n)
    for script in scripts:
        i = 0
        j = 0
        inDiag = 0
        for s in script:
            i_next = i
            j_next = j
            if (s == 'M') or (s == 'S') or (s == 'D'):
                i_next += 1
            if (s == 'M') or (s == 'S') or (s == 'I'):
                j_next += 1
            # enter diagonal
            if ( (j_next - i_next) == diag ) and ((j - i) != diag):
                enterStats[i_next] += 1
                inDiag += 1
            # exit diagonal
            if ((j_next - i_next) != diag) and ((j - i) == diag):
                exitStats[i] += 1
                onDiagStats[inDiag] += 1
                inDiag = 0
            i = i_next
            j = j_next
    return (enterStats, exitStats, onDiagStats)

def diagonalStatAverageLengths(n, scripts, diagonals):
    print("Hello!!")
    for diag in diagonals:
        ent, ex, lens = calculateDiagonalStats(n, scripts, diag)
        normLengths = [t / sum(lens) for t in lens]
        indexes = [t for t in range(1, len(lens)+1)]
        avg = sum([(i+1) * lens[i] / sum(lens) for i in range(len(lens))])
        print ("{0}\t{1}".format(diag, avg))
                                 #np.average(normLengths, weights=indexes)))
    

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
    mDist, sDist, dDist, iDist = operationsDistribution(n+1, allScripts)

    # WARNING: Using python3 the savetxt with zip fails    
    # pack into one file distributions vector for all operations along with the header
    if sys.version_info.major <= 2:
        np.savetxt(args.distFile ,zip(np.arange(n+1), mDist, sDist, dDist, iDist), delimiter="\t")
    else:
        print("{0} {1} file not produced (Python{2} issue)".format(yellowWarnMsg, args.distFile, sys.version_info.major))
    np.savetxt(args.matFile, freqMat, delimiter="\n")
    np.savetxt(args.maxFile, maxDistr, delimiter="\n")
    # cross distribution requested
    flow_r = None
    flor_c = None
    r = -1
    c = -1
    if args.matrixRow:
        r = int(args.matrixRow)
        flow_r = calculateHorizontalFlow(n+1, allScripts, r)
        for i in range(n+1):            
            print( "{0}\t{1}".format( i, int(flow_r[i]) ) )
    if args.matrixColumn:
        c = int(args.matrixColumn)
        flow_c = calculateVerticalFlow(n+1, allScripts, c)
        for j in range(n+1):
            print( "{0}\t{1}".format( j, int(flow_c[j]) ) )

    if (args.stretchStat is not None):
        mStretch, sStretch, dStretch, iStretch = computeStretchesDistribution(n+1, allScripts)

    if(args.diagonalStats is not None):
        diagEnter, diagExit, diagLength = calculateDiagonalStats(n+1, allScripts)
        diagonalStatAverageLengths(n+1, allScripts, [-2,-1,0,1,2])
            
    # The archive options puts all files in binary numpy format into
    # the archiviation directory. This is useful to load afterwards
    # with the ndarray.fromfile function (useful e.g., for plotting).
    # The name of the files are derived from the corresponding plain
    # text (given as parameter or defualted).
    if args.archive:
        archivePath = path.abspath(args.archive)
        existPath = path.exists(archivePath)
        # indicated path is a file, can not proceed
        if existPath and (not path.isdir(archivePath)):
            print( "{0} {1} (option -a) is not a directory!".format(redErrorMsg, args.archive))
            sys.exit(1)
        if not existPath:
            os.mkdir(archivePath)
        np.save(os.path.join(archivePath, args.matFile), freqMat )
        # make one file for each type of operation to make loading easire
        np.save(os.path.join(archivePath, "M_" + args.distFile), mDist)
        np.save(os.path.join(archivePath, "S_" + args.distFile), sDist)
        np.save(os.path.join(archivePath, "D_" + args.distFile), dDist)
        np.save(os.path.join(archivePath, "I_" + args.distFile), iDist)
        np.save(os.path.join(archivePath, args.maxFile), maxDistr)
        # save also crossing distributions if calculated
        if args.matrixRow:
            np.save(os.path.join(archivePath, "HCross_" + str(r)), flow_r)
        if args.matrixColumn:
            np.save(os.path.join(archivePath, "VCross_" + str(c)), flow_c )

        # save stretches
        if (args.stretchStat is not None):
            np.save(os.path.join(archivePath, "M_" + args.stretchStat), mStretch)
            np.save(os.path.join(archivePath, "S_" + args.stretchStat), sStretch)
            np.save(os.path.join(archivePath, "D_" + args.stretchStat), dStretch)
            np.save(os.path.join(archivePath, "I_" + args.stretchStat), iStretch)
    
        if (args.diagonalStats is not None):
            np.save(os.path.join(archivePath, "Enter_" + args.diagonalStats), diagEnter)
            np.save(os.path.join(archivePath, "Exit_" + args.diagonalStats), diagExit)
            np.save(os.path.join(archivePath, "Len_" + args.diagonalStats), diagLength)

