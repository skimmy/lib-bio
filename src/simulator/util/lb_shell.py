"""This module contains all functions needed to evaluate a
shell-based lower bound to the limit constant for the expected edit
distance 

Todo:
   * functions for calculating rstar (normalized and not normalized)
   * extensive tests
   * functions for our bounds: insertion distance, simple bound and
     canonical path based bound"""

import numpy as np
import scipy.special as spec

def bound_with_shells(n, rstar, sigma, shell_ub, norm=True):
    """This function returns lower bound to \alpha_n given upper bounds to
    shells (see eq. 7 of the report).

    Parameters
    ----------
    n : strings length
    rstar : maximum radius
    sigma : alphabet size
    shell_ub : a function accepting three integers r, n, sigma and returns an upper
       bound to the size of the radius r shell
    norm : indicates whether or not the upper bound is normalized in [0,1],
       if False, the function normilizes it dividing by sigma^n.
       It is important to use normalized function when dealing with big n
       otherwise the computation of the bound may overflow.

    Returns
    -------
    A lower bound to e(n)/n
      

    """
    s = 0
    norm_fact = 1
    if(not norm):
        norm_fact = np.power(sigma,n) 
    for r in range(rstar):
        for i in range(r):
            s += shell_ub(i, n, sigma) / norm_fact
    return (rstar - s)/n

def bound_with_volumes(n, rstar, sigma, vol_ub, norm=True):
    """This function returns lower bound to \alpha_n given upper bounds to
    volumes (see eq. 8 of the report)
    Parameters
    ----------
    n : strings length
    rstar : maximum radius
    sigma : alphabet size
    vol_ub : a function accepting three integers r, n, sigma and returns an upper
       bound to the size of the radius r volume
    norm : indicates whether or not the upper bound is normalized in [0,1],
       if False, the function normilizes it dividing by sigma^n.
       It is important to use normalized function when dealing with big n
       otherwise the computation of the bound may overflow.

    Returns
    -------
    A lower bound to e(n)/n
      

    """
    s = 0
    norm_fact = 1
    if(not norm):
        norm_fact = np.power(sigma,n) 
    for r in range(rstar):
        s += vol_ub(r, n, sigma) / norm_fact
    return (rstar - s)/n

def find_rsat_volume_n(n, sigma, v_n, step=1):
    rs = 0
    v = v_n(rs, n, sigma)
    while(v<1):
        rs += step
        v = v_n(rs, n, sigma)
    return rs-step

# Insertion distance bounds
def vr_ub_n_2(r, n, sigma):
    return spec.binom(n,r)**2 / np.power(sigma, n-r)

def vr_ub_n(r, n, sigma):
    p = 0
    log2_sigma = np.log2(sigma)
    for i in range(1,n-r+1):
        p += 2*np.log2(1+r/i) -  log2_sigma
    return np.power(2,p)

# Simple bound
def wr_ub_n_r(r, n, sigma):
    s = 0
    dmax = int(np.floor(r/2))
    for d in range(dmax+1):
        s += (spec.binom(n,d)**2)*spec.binom(n-d,r-2*d)*np.power(4/9,d)
    return s*np.power(3,r)/np.power(sigma,n)

def util_sum_log2(low, high):
    s = 0
    for i in range(low,high+1):
        s += np.log2(i) # could use pre-computed table here
    return s

def wr_ub_n_r_2(r, n, sigma):
    dmax = int(np.floor(r/2))
    s = 0
    for d in range(dmax+1):
        C = 1 + np.log2(3.0)*(r-2*d) + 2*(d-n)
        L1 = 0
        for i in range(1,d+1):
            L1 += np.log2((n-i+1)/i)
        L2 = 0
        for i in range(1,r-2*d+1):
            L2 += np.log2((n-d-i+1)/i)
        L = 2*L1 + L2
        s+= np.power(2,C + L)
    return s

def wr_ub_n_2(r, n, sigma):
    if (r == 0):
        return np.power(sigma,-float(n))
    wr = wr_ub_n_r_2(r, n, sigma)
    wr_1 = wr_ub_n_r_2(r-1, n, sigma)
    return wr + wr_1

def wr_ub_n(r, n, sigma):
    # notice that the bound is W_r + W_{r-1}
    if (r == 0):
        return 1/np.power(sigma,n)
    wr = wr_ub_n_r(r, n, sigma)
    wr_1 = wr_ub_n_r(r-1, n, sigma)
    return wr+wr_1

def evaluate_insertion_distance_bound():
    sigma = 4
    N = 10000
    for n in range(401,N,239):
        rs = find_rsat_volume_n(n, sigma, vr_ub_n, step=7)
        v = bound_with_volumes(n, rs, sigma, vr_ub_n)
        print("{0}\t{1}\t{2:.4f}".format(n, rs, v, step=1))

def evaluate_simple_bound():
    sigma = 4
    N = 10000
    for n in range(10,N,100):
        rs = find_rsat_volume_n(n,sigma, wr_ub_n_2, step=1)
        bound = 0 #bound_with_volumes(n, rs, sigma, wr_ub_n)
        bound2 = bound_with_volumes(n, rs, sigma, wr_ub_n_2)
        print("{0}\t{1}\t{2}".format(n, bound, bound2))

if (__name__ == "__main__"):
    #evaluate_insertion_distance_bound()
    evaluate_simple_bound()
