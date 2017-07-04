import numpy as np
from scipy.stats import binom

def expectation(pmf):
    return sum(np.multiply(pmf, np.arange(len(pmf))))

def binomial_distribution_vectors(n, p, v_size=None):
    pmf = np.array([binom.pmf(i, n, p) for i in range(n+1)])
    cdf = np.array([binom.cdf(i, n, p) for i in range(n+1)])
    if (v_size != None) and (v_size > n+1):
        # padding to the right to reach v_size length
        pmf = np.append(pmf, np.zeros(v_size-n))
        cdf = np.append(cdf, np.ones(v_size-n))
    return pmf,cdf

def minimum_distribution_vectors(X, Y, n):
    pmf = np.zeros(n+1)
    for i in range(n+1):
        p = 0
        # X,Y equals
        p += X[i]*Y[i]
        # min is X
        p += X[i] * sum(Y[i+1:])
        # min is Y
        p += Y[i] * sum(X[i+1:])
        pmf[i] = p
    cdf = np.array([sum(pmf[0:i+1]) for i in range(n+1)])
    return pmf, cdf
    


if (__name__ == "__main__"):
    p = 0.75
    for n in range(4,65):
        # shift distribution support: 2,3,...,n+1 (padding to the left)
        binom_t_1, c_t_1 = binomial_distribution_vectors(n-1, p, n)
        shift_pmf = np.append(np.zeros(2), binom_t_1)
        # non-shift distribution support 0,1,...,n (padding to the right)
        p_mass, cum_dist = binomial_distribution_vectors(n, p, n+2)

        p_min, cum_min = minimum_distribution_vectors(shift_pmf, shift_pmf, n+1)        
        p_min, cum_min = minimum_distribution_vectors(p_min, p_min, n)
        e_val = expectation(p_min)
        print("{0}\t{1}\t{2}\t{3}\t{4}".format(n, e_val, p*n, p*n-e_val, e_val / n ))
