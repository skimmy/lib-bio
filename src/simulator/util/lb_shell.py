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

def bbb(x, n, sigma=4):
    return np.power(sigma, x)

def bbb_n(x, n, sigma=4):
    return np.power(sigma,x) / np.power(sigma,n)

def vr_ub_n_2(r, n, sigma):
    return spec.binom(n,r)**2 / np.power(sigma, n-r)

def vr_ub_n(r, n, sigma):
    p = 0
    log2_sigma = np.log2(sigma)
    for i in range(1,n-r+1):
        p += 2*np.log2(1+r/i) -  log2_sigma
    return np.power(2,p)

if (__name__ == "__main__"):
    # n = 10
    # rs = 3
    # sigma = 4
    # a = bound_with_shells(n, rs, sigma, bbb, False)
    # an = bound_with_shells(n, rs, sigma, bbb_n)
    # print(a)
    # print(an)

    # b = bound_with_volumes(n, rs, sigma, bbb, False)
    # bn = bound_with_volumes(n, rs, sigma, bbb_n)
    # print(b)
    # print(bn)
    sigma = 4
    N = 10000
    for n in range(401,N,239):
        rs = find_rsat_volume_n(n, sigma, vr_ub_n, step=7)
        v = bound_with_volumes(n, rs, sigma, vr_ub_n)
#        v2 = bound_with_volumes(n, rs, sigma, vr_ub_n_2)
        print("{0}\t{1}\t{2:.4f}".format(n, rs, v, step=1))
