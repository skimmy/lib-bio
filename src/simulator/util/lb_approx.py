import sys
import math
import scipy.special as spec

def s_r_UB_q_approx(n, r, nD, q, Sigma):
    delta_term = ( spec.binom(n-r+nD,q) + spec.binom(n-r+nD,q-1) ) 
    return delta_term * spec.binom(q , int(math.floor(q/2)))*spec.binom(2*nD-1, q-2)

def s_r_UB_q(n, r, nD, q, Sigma):
    return s_r_UB_q_approx(n, r, nD, q, Sigma)

def s_r_UB_nD(n, r, nD, Sigma):
    sigma_factor = (pow(Sigma-1, nD) + pow(Sigma, nD)) / pow(Sigma-1, 2*nD)
    q_sum = 0
    q_max = min(2*nD, n - r + nD + 1) # bar{q} = min{2nD, nM+1}
    for q in range(2,q_max+1):
        q_sum += s_r_UB_q(n, r, nD, q, Sigma)
    return sigma_factor* spec.binom(n-nD, r-2*nD) * (n / (n-nD)) * q_sum

def s_r_UB(n, r, Sigma):
    nD_sum = 0
    if (r == 0):
        nD_sum = spec.binom(n,r)
    else:
        nD_max = int(math.floor(r/2))
        for nD in range(nD_max+1):
            nD_sum += s_r_UB_nD(n, r, nD, Sigma)
    return nD_sum*pow(Sigma-1, r)

def s_UB(n, Sigma):
    r_max = n
    shell_bounds = []
    for r in range(r_max+1):
        shell_bounds.append(s_r_UB(n, r, Sigma))
    return shell_bounds

def bound(n, Sigma=4): # OK
    hulls_bound = [(i,0) for i in range(n+1)]    
    remaining_strings = pow(Sigma,n) - 1
    lb = 0
    r = 1    
    while (remaining_strings > 0 and r <= n):
        lb += remaining_strings
        hull_count = s_r_UB(n,r, Sigma)             
        remaining_strings = remaining_strings - hull_count
        hulls_bound[r] = (r,lb)
        r += 1
    
    return (lb, hulls_bound, r-1)

if (__name__ == "__main__"):
    n_max = int(sys.argv[1])
    n_step = int(sys.argv[2])
    Sigma=4
    print("n_max: {0}\nn_step: {1}".format(n_max, n_step))
    ns = list(range(50,n_max+1,n_step))
    for n in ns:
        UB = 0
        shell_bounds = s_UB(n, Sigma)        
        (lb, hulls_bound, r_sat) = bound(n)
        print("{0}\t{1}".format(n,lb / (n*pow(Sigma,n))))
        
        
