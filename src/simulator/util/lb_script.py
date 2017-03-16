import math
import scipy.special as spec

'''Convenienve function (if ever needed to improve) to compute
  binom(a,b) * binom(c,d) * Sigma^d'''
def Fd(a,b,c,d,Sigma):
    return spec.binom(a,b) * spec.binom(c,d) * pow(Sigma-1,d)

def compute_fd00(n, r, D, q, Sigma):
    # the case q=M+1=n-r+d+1 is not possible
    if (q >= n-r+D+1):
        return 0
    return Fd(n-D+1, q, n-D-q, r-2*D, Sigma)

def compute_fd01(n, r, D, q, Sigma):
    return 0

def compute_fd10(n, r, D, q, Sigma):
    # the case q=M+1=n-r+d+1 is not possible
    if (q >= n-r+D+1):
        return 0
    return 0

def compute_fd11(n, r, D, q, Sigma):    
    return 0

'''This function computes the number of scripts that can be generated
through a sufficient path for a given number of non diagonal path q
and number of Deletion/Insertion D'''
def compute_sufficient_path_D(n, r, D, q, Sigma):
    # q=0 there are only diagonal arcs
    if (q==0):
        return ( spec.binom(n,r) * pow(Sigma-1,r) )
    fd00 = compute_fd00(n, r, D, q, Sigma)
   # print("n={0} r={1} D={2} q={3}\t\tfd00={4}".format(n,r,D,q,fd00))
    fd01 = compute_fd01(n, r, D, q, Sigma)
    fd10 = compute_fd10(n, r, D, q, Sigma)
    fd11 = compute_fd11(n, r, D, q, Sigma)
    return fd00


'''Compute the volume V_r under the constraint that exactly D
Deletion and Insertion operation are used'''
def compute_volume_r_D(n, r, D, Sigma):
    # Number of matches (n-r-
    M = n - r + D
    # q = 0,2,3,...,min(2D,M+1)
    q_min = 2
    q_max = min(2*D, M+1)
    q_range = (list(range(q_min,q_max+1)))
    q_range.insert(0,0)
    Q_sum = 0
    for q in q_range:
        Q_sum += compute_sufficient_path_D(n, r, D, q, Sigma)
    return Q_sum


'''For a given number of errors  0 <= r <= n computes an upper bound
V_r to the number of strings at distance <= r from any center A'''
def compute_volume_r_UB(n, r, Sigma):
    D_max = r >> 1
    Vr = 0
    for D in range(D_max+1):
        Vr += compute_volume_r_D(n, r, D, Sigma)
#    return pow(Sigma,n)
    return Vr

'''Compute the lower bound on constant alpha for given n'''
def compute_lower_bound(n, Sigma=4):
    norm_const = float(pow(Sigma, n) * n)
    # The lower bound is obtained via an upper bound on the volume
    volumes_sum = 0
    # any value in [1,n-1] contributes to the bound (see LaTeX)
    for r in range(1,n):
        vr = compute_volume_r_UB(n, r, Sigma)
        volumes_sum += vr

        # consider moving normalization inside the above for loop
    return 1 - (float(volumes_sum) / norm_const)

if __name__ == "__main__":
    n_max = 32
    for n in range(1,n_max+1):
        print("{0}\t{1}".format(n,compute_lower_bound(n)))
