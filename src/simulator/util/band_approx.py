#!/usr/bin/python

import argparse
import numpy as np

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="File with results of simulation")
    parser.add_argument("n", help="length of the strings")
    parser.add_argument("-o", dest="outfile", help="Name of the output file", default="out.txt")
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

'''Returns a dictionary with a numpy vector for each of the entries
of Ts. Such vector will contain in position 'i' the number of times
a difference 'i' was observed between exact and approximated version'''
def computeDistributions(data, Ts, n):    
    dists = {}
    for T in Ts:
        dists[T] = np.zeros(n)
    for d in data:
        exact = sum(d[0])
        for T in Ts:
            diff = sum(d[T]) - exact
            dists[T][diff] += 1
    return dists

def saveDistributionsToFile(dists, filename, n):
    Ts = [T for T in dists if T > 0]
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
    
    fh.close()
