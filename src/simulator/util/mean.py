#!/usr/bin/python

# Compute the mean of values from the input file 

import argparse

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="File with value to mean")
    # parser.add_argument("-o", "--optional", dest='o', help="Optional flag", action='store_true')
    # parser.add_argument("-d", "--default", help="SWith default", default="Hello")
    return parser.parse_args()

def computeMean(l):
    return (sum(l) / float(len(l)))
    

if __name__ == "__main__":
    args = parseArguments()
    fh = open(args.input_file)
    values_list = []
    for line in fh:
        values_list.append(float(line))
    print(computeMean(values_list))
    fh.close()
