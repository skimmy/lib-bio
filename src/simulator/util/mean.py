#!/usr/bin/python

# Compute the mean of values from the input file 

import argparse

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="File with value to mean")
    parser.add_argument("-V", "--variancel", dest='V', help="Computes variance", action='store_true')
    # parser.add_argument("-d", "--default", help="SWith default", default="Hello")
    return parser.parse_args()

def computeMean(l):
    return (sum(l) / float(len(l)))

def computeVariance(l):
    n = len(l)
    v = [0]*n
    mean = computeMean(l)
    for i in range(n):
        v[i] = (l[i] - mean) * (l[i] - mean)
    return computeMean(v)
    

if __name__ == "__main__":
    args = parseArguments()
    fh = open(args.input_file)
    values_list = []
    for line in fh:
        values_list.append(float(line))
    if args.V:
        print(computeVariance(values_list))
    else:
        print(computeMean(values_list))
    
    fh.close()
