#!/usr/bin/python

import argparse
import numpy as np
import math

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="File with results of simulation")
    parser.add_argument("n", help="length of the strings")
    parser.add_argument("-o", dest="outfile", help="Name of the output file", default="out.txt")
    parser.add_argument("-P","--partition", dest="partition",help="Name of the file for operation distribution", default=None)
    # parser.add_argument("opt", help="Required option")
    # parser.add_argument("-o", "--optional", dest='o', help="Optional flag", action='store_true')
    # parser.add_argument("-d", "--default", help="SWith default", default="Hello")
    return parser.parse_args()

def printData(data, Ts):    
    print("\t\t".join([str(x) for x in Ts]))
    for d in data:
        s = ""
        for i in range(len(Ts)):
            s += " ".join([str(x) for x in d[Ts[i]] ])
            s += "\t"
        print s

def printDistances(data, Ts):
    for d in data:
        print("\t".join([ str(sum( d[Ts[i]] )) for i in range(len(Ts))]))

def printStatistics(params, Ts, n):
    vals = [T for T in Ts if T < n]
    list.sort(vals)
    for v in vals:
        print("{0}\t{1}\t{2}\t{3}".format(v, params[v][0], math.sqrt(params[v][1]), params[v][2]))

def printOperationsStatistics(part, Ts):
    vals = [T for T in Ts]
    list.sort(vals)
    for T in vals:
        meanS = 0
        meanD = 0
        meanI = 0
        for x in range(len(part[T]['S'])):
            meanS += x * part[T]['S'][x]
            meanI += x * part[T]['D'][x]
            meanD += x * part[T]['I'][x]
        meanS /= sum(part[T]['S'])
        meanD /= sum(part[T]['D'])
        meanI /= sum(part[T]['I'])
        print("{0}\t{1}\t{2}\t{3}".format(T,meanS, meanD, meanI))
    

def saveDistributionsToFile(dists, filename, n):
    Ts = [T for T in dists if T < n]
    list.sort(Ts)
    print Ts
    fh = open(filename, "w")
    for i in range(n):
        s = str(i) + "\t"
        for T in Ts:
            s += str(int(dists[T][i])) + "\t"
        s += "\n"
        fh.write(s)
    fh.close()

def saveOperationsPartitionToFile(part, Ts, fileName):
    fh = open(fileName, "w")
    for T in Ts:
        pass
    fh.close()


'''Returns a dictionary containing an entry for each of the value of Ts. Each of
such entry contains a dictionary with 3 entries indexed by 'S', 'D' and 'I' and
containing (as value) a numpy vector V of size n such that V[i] contains the number
of times of substitution, deletion and insertion (respectively) observed.

Calling P the returned dictionary the element
  P[T]['S'][i]
is the number of times 'i' substitution ('S') where observed in all the scripts
obtained using the approximation value T
'''
def computeOperationsPartition(data, Ts, n):
    part = {}
    for T in Ts:
        part[T] = {}
        part[T]['S'] = np.zeros(n)
        part[T]['D'] = np.zeros(n)
        part[T]['I'] = np.zeros(n)
        for d in data:
            part[T]['S'][(d[T][0])] += 1
            part[T]['D'][(d[T][1])] += 1
            part[T]['I'][(d[T][2])] += 1
    return part


        
'''Returns a dictionary with a numpy vector for each of the entries
of Ts. Such vector will contain in position 'i' the number of times
a difference 'i' was observed between exact and approximated version'''
def computeDistributions(data, Ts, n):    
    dists = {}
    for T in Ts:
        dists[T] = np.zeros(n)
    for d in data:
        exact = sum(d[n])
        for T in Ts:
            diff = sum(d[T]) - exact
            dists[T][diff] += 1
    return dists

def computeStatistics(dists, Ts, n):
    params = {}
    for T in Ts:
        freqs = dists[T]
        k = sum(freqs)
        expectation = 0
        variance = 0
        maxDiff = 0
        for i in range(len(freqs)):
            pi = freqs[i] / k
            expectation += pi * i
            variance += pi * i * i
            if pi > 0:
                maxDiff = i
        params[T] = (expectation, (variance - expectation*expectation), maxDiff)
    return params

if __name__ == "__main__":
    args = parseArguments()
    n = int(args.n)
    fh = open(args.file, "r")
    data = []
    apprValues = fh.readline().strip().split("\t")
    Ts = [int(x) for x in apprValues]
    for line in fh:
        entries = line.strip().split("\t")
        entry = {}
        for i in range(len(entries)):
            info = [int(x) for x in entries[i].split(" " )]
            entry[Ts[i]] = info
        data.append(entry)
    dists = computeDistributions(data, Ts, n)
    saveDistributionsToFile(dists, args.outfile, n)
    params = computeStatistics(dists, Ts, n)
    printStatistics(params,Ts, n)    
    if args.partition != None:
        print(Ts)
        part = computeOperationsPartition(data, Ts, n)
        saveOperationsPartitionToFile(part, Ts, args.partition)
        printOperationsStatistics(part, Ts)
    fh.close()
