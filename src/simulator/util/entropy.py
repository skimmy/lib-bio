#!/usr/bin/python

import sys
import math
import heapq
import scipy.special as spe

# calculates the product of elements in v
#    p = \prod_{f \in v}{f}
def product(v):
    p = 1.0;
    for f in v:
        p = p*f
    return p

# Inverts all the signs of elements of v
#    v[i] = -v[i]
def listSignInvert(v):
    for i in range(len(v)):
        v[i] = -v[i]

# This function calculates multinomial coefficient of
# ks = [k1, k2, ..., kn] by keeping the factorized versione
# of numerator and denumerator. This allows the calculation
# not to overflow as long as the final results doesn't
#
# This procedure makes use of heap structure implemented in the
# python module 'heapq'.
# Notice that it returns two lists one containing factors of the
# numerator and the other the factors of denumerator. This latter
# list, upon success, is always empty. Also notice that factors
# may be floating numbers and, therefore, the product may be a
# float which (for all practical cases) can be rounded to an integer
# representing the wanted result.
def multinomial(ks):
    n = len(ks)
    if n <= 1:        
        return ([1.0],[]);
    k = sum(ks)
    k1 = ks[0]
    h_num = []
    for i in range(k1+1,k+1):
        heapq.heappush(h_num, -i)
    
    h_den = []
    for j in range(1,n):
        ki = ks[j]
        for i in range(1,ki+1):
            heapq.heappush(h_den,-i)
    h_den = sorted(h_den)
    while h_den:
        max_num = -heapq.heappop(h_num)
        max_den = -heapq.heappop(h_den)
        new_fact = float(max_num) / float(max_den)
        heapq.heappush(h_num,-new_fact)
    listSignInvert(h_num)
    listSignInvert(h_den)
    return (h_num,h_den)



# Returns a list containing all the partitions of integer 'k'
# into 'n' integers (possibly with repetitions)
def recPart(k,n,kmax=-1):
    if (kmax==-1):
        kmax = k+1        
    if (k==0):
        return [[0]*n]
    if (n==1):
        if (k > kmax):
            return []
        return [[] +[k]]
    kbar = int(math.ceil(k/n))
    P = []
    kmax = min(k,kmax)
    for i in range(kmax,kbar-1,-1):       
        Q = recPart(k-i,n-1,i)
        for q in Q:            
            P.append([i] + q)
    return P

def multiplicities(ks):
    mus = []
    current = ks[0]
    mu  = 1
    for i in range(1,len(ks)):
        if current == ks[i]:
            mu += 1
        else:
            mus.append(mu)
            mu = 1
            current = ks[i]
    mus.append(mu)
    return mus

# Tests the partition algorithm by summing all 4^k observations
def testPartition(P):
    s = 0
    
    for p in P:
        (c1,tmp1) = (multinomial(p))    
        (c2,tmp2) = multinomial(multiplicities(p))        
        s += product(c1) * product(c2)
    return s

# Main procedure
# Usage
#    npart k [n]
#
#      k    integer to be partitioned
#      n    number of parts (optional, default value: 4)
if __name__ == "__main__":
    n = 4
    k = int(sys.argv[1])
    if len(sys.argv) > 2:
        n = int(sys.argv[2])
    P = recPart(k,n)
    print("There are {0} {1}-partitions of {2}:").format(len(P),n,k)
    print(testPartition(P))
    print(math.pow(4,k))
    
