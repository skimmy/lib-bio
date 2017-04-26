import random
import re
import sys

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
    n = 12
    r = 5
    D = 2
    k = 5
    if (len(sys.argv) > 1):
        k = int(sys.argv[1])
    do_print = (sys.argv.count("-print"))

    # construct Regex for canonical path
    pattern  = r'(S*M*((D+|I+)M)*)*(D*|I*)$'
    matched = 0
    for i in range(k):
        path = "".join(generate_random_string(n,r,D))
        if (re.match(pattern, path)):
            matched += 1
            if (do_print):
                print("+ {0}".format(path))
        else:
            if(do_print):
                print("- {0}".format(path))

    print("Mathced {0}/{1} Ratio: {2}".format(matched, k, matched/k))
