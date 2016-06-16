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
def ListSignInvert(v):
    for i in range(len(v)):
        v[i] = -v[i]

def ListRemoveZeros(v):
    return [x for x in v if x != 0]

def ListRemoveOnes(v):
    return [x for x in v if x != 1]

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
    ListSignInvert(h_num)
    ListSignInvert(h_den)
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

########## MULTINOMIAL FACTORS PROCEDURES ##########
def MultinomialFactors(ks):
    q = [] + ListRemoveOnes(ListRemoveZeros(ks))
    ListSignInvert(q)
    heapq.heapify(q)
    return [sum(ks), q]

def Next(Q):
    if Q[0] <= 1:
        return 1
    d = 1
    if (Q[1]):
        d = -heapq.heappop(Q[1])
        if d-1 > 1:
            heapq.heappush(Q[1],-(d-1))
    Q[0] = Q[0] - 1    
    return float(Q[0] + 1) / float(d)
            
def IsEmpty(Q):
    return (Q[0] <= 1 and (not Q[1]))

def SafeProbabilityComputing(ks, eps):
    k = sum(ks)
    S = 0
    for i in range(4):
        Q = MultinomialFactors(ks)
        P = 1.0
        ki = ks[i]
        # go trhough (1-eps)
        for i in range(ki):
            P *= 1.0-eps
            while(P < 1) and (not IsEmpty(Q)):
                P *= Next(Q)
        # go through eps/3
        for i in range(k-ki):
            P *= eps/3.0
            while(P < 1) and (not IsEmpty(Q)):
                P *= Next(Q)
        S += P
    mu = 0.25 * product(multinomial(multiplicities(ks))[0] )
#    print("{0} {1} {2} {3}").format(S,mu, product(multiplicities(ks)), ks)
    return S *mu


########## ENTROPY CALCULATION ##########
def SimplePartitionEntropy(ks, eps):
    H = 0    
    for s in range(len(ks)):
        etas = 0;
        for t in range(len(ks)):
            if (t != s):
                etas += math.pow(1-eps,ks[t]-ks[s]) * math.pow(eps/3, ks[s]-ks[t])
        H += 1 / (1 + etas) * (math.log(1 + etas) / math.log(2) )
    return H
   


def TestMultinomialFactorsCalculator(ks):
    Q = MultinomialFactors(ks)   
    print(Q)
    p = 1
    while (not IsEmpty(Q)):
        x = Next(Q)
        p *= x
        print("{0} - {1}").format(x,p)
    mn = product(multinomial(ks)[0])
    return (int(mn) == int(p))

def TestSimpleEntropy(ks, eps):
    print(SimplePartitionEntropy(ks, eps))

       
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
    epsilon = 0.21
    n = 4
    k = int(sys.argv[1])
    if len(sys.argv) > 2:
        n = int(sys.argv[2])
    prob = 0
    entropy = 0
    P = recPart(k,n)
    #    print("There are {0} {1}-partitions of {2}:").format(len(P),n,k)
    #    p = [10,1,0,0]
    #    p = [1,1,1,1]
    #    print(TestMultinomialFactorsCalculator(p))
    #    TestSimpleEntropy([2,1,1,0],0.21)
    for p in P:
        pr = float(SafeProbabilityComputing(p,epsilon))
        prob += pr
        entropy += pr * SimplePartitionEntropy(p, epsilon)
    print("{0}").format(entropy)

