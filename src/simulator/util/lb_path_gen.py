import random
import re
import sys
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser()
    # parser.add_argument("opt", help="Required option")
    parser.add_argument("-p", "--pipeline", help="When set only shows result ratio", action='store_true')
    parser.add_argument("-n", "--length", help="Length of the strings n", default=16)
    parser.add_argument("-r", "--errors", help="Total number of errors r", default=8)
    parser.add_argument("-D", "--deletions", help="Number of insertions and deltions D", default=4)
    parser.add_argument("-k", "--samples", help="Number of samples k", default=32)

    return parser.parse_args()


def generate_random_string(n, r, D):
    script = ['M']*(n+D)
    l = [x for x in range(n+D)]
    # add Indels
    for i in range(D):
        r = random.choice(l)
        l.remove(r)
        script[r] = 'D'
        r = random.choice(l)
        l.remove(r)
        script[r] = 'I'
    for i in range(r-2*D):
        r = random.choice(l)
        l.remove(r)
        script[r] = 'S'
    return script


if (__name__ == "__main__"):
    args = parse_arguments()
    n = int(args.length)
    r = int(args.errors)
    D = int(args.deletions)
    k = int(args.samples)
    do_print = (sys.argv.count("-print"))

    # construct Regex for canonical path
    pattern  = r'(S*M*((D+|I+)M)*)*(D*|I*)$'
    matched = 0
    for i in range(k):
        path = "".join(generate_random_string(n,r,D))
        if (re.match(pattern, path)):
            matched += 1

    if (args.pipeline):
        print(matched/k)
    else:
        print("Mathced {0}/{1} Ratio: {2}".format(matched, k, matched/k))
