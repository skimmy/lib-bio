import math
import scipy.special as spec

##############################################################################
#                                                                            #
#                         CONVENIENCE FUNCTIONS                              #
#                                                                            #
##############################################################################

'''Convenienve function (if ever needed to improve) to compute
  binom(a,b) * binom(c,d) * Sigma^d'''
def Fd(a,b,c,d,e):
    return spec.binom(a,b) * spec.binom(c,d) * pow(e,d)

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

# Case q=0 (i.e., main diagonal path)
def computeQ0(n, r, Sigma):
    return spec.binom(n, r) * pow(Sigma-1, r)

def compute_fd(n, r, D, q, delta_s, delta_e, Sigma):    
    t = q + 1 - (delta_s + delta_e) # t = q + delta
    return Fd(n-D-1, t-1, n-D-q+delta_e, r - 2*D, Sigma-1)
    
def compute_fn(n, r, D, q, delta_s, delta_e, Sigma):
    count = 0
    # no diagonal segment at the end
    if (delta_e == 0):
        d_max = min(q-1,D)
        for d in range(1,d_max+1):
            count += (spec.binom(q,d) * spec.binom(D-1, d-1)
                      * spec.binom(D-1,q-d-1))
        count *= pow(Sigma-1, D)
    # diagonal at the end
    else:
        # ending del
        d_bar = min(q-1,d)
        for d in range(d_bar):
            count += (spec.binom(q-1,d) * spec.binom(D-1,d)
                      * spec.binom(D-1, q-d-1))
        count *= pow(Sigma-1, D)
        
        # ending ins
        i_bar = min(q-1, D)
        for i in range(i_bar):
            part_count = 0
            for j in range(1,D-i):
                part_count += (spec.binom(D-j-1,i-1) * spec.binom(D-1,q-i) *
                               pow(Sigma,j) * pow(Sigma-1,D-j))
            count += part_count * spec.binom(q-1,i)
            
    return count

##############################################################################
#                                                                            #
#                          OLD BOUND FUNCTIONS                               #
#                                                                            #
##############################################################################

'''This function computes the number of scripts that can be generated
through a sufficient path for a given number of non diagonal path q
and number of Deletion/Insertion D'''
def compute_sufficient_path_D(n, r, D, q, Sigma):
    fd00 = compute_fd00(n, r, D, q, Sigma)
    fd01 = compute_fd01(n, r, D, q, Sigma)
    fd10 = compute_fd10(n, r, D, q, Sigma)
    fd11 = compute_fd11(n, r, D, q, Sigma)
    return fd00


'''Compute the volume V_r under the constraint that exactly D
Deletion and Insertion operation are used'''
def compute_volume_r_D(n, r, D, Sigma):
    # Number of matches (n-r+D)
    M = n - r + D
    q_max = min(2*D, M+1)
    Q_sum = 0
    # q = 2,3,...,min(2D,M+1)
    for q in range(2,q_max):
        Q_sum += compute_sufficient_path_D(n, r, D, q, Sigma)
    return Q_sum


'''For a given number of errors  0 <= r <= n computes an upper bound
V_r to the number of strings at distance <= r from any center A'''
def compute_volume_r_UB(n, r, Sigma):
    D_max = r >> 1
    Vr = computeQ0(n, r, Sigma)
    for D in range(D_max+1):
        Vr += compute_volume_r_D(n, r, D, Sigma)
    return Vr 

'''Computes the lower bound on constant alpha for given n'''
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

##############################################################################
#                                                                            #
#                          NEW BOUND FUNCTIONS                               #
#                                                                            #
##############################################################################

def count_canonical_annotated_path_r_D_q(n, r, D, q, Sigma):
    count = 0
    # sum of fd*fn
    for delta_s in range(1):
        for delta_e in range(1):
            fd = compute_fd(n, r, D, q, delta_s, delta_e, Sigma)
            fn = compute_fn(n, r, D, q, delta_s, delta_e, Sigma)
            count += fd*fn
    return count
            
def count_canonical_annotated_path_r_D(n, r, D, Sigma):
    # TODO ...
    count = 0
    q_max = min(2*D, n - r + D) ####### CHECK!!!!!!!
    for q in range(2,q_max+1):
        count += count_canonical_annotated_path_r_D_q(n, r, D, q, Sigma)
    return count

'''Counts the number of canonical annotated paths of cost r'''
def count_canonical_annotated_path_r(n, r, Sigma):
    count = computeQ0(n, r, Sigma) # Q0 + ...
    Dmax = int(math.floor(r / 2.0)) # !!!!!! CHECK !!!!!!
    for D in range(1,Dmax+1):
        count += count_canonical_annotated_path_r_D(n, r, D, Sigma)
    return count;

def lower_bound(n, Sigma=4):
    lb = 1
    for r in range(1,n+1):
        lb += r * count_canonical_annotated_path_r(n, r, Sigma)
    return lb / pow(Sigma, 2*n)

##############################################################################
#                                                                            #
#                                 MAIN                                       #
#                                                                            #
##############################################################################

if __name__ == "__main__":
    n_max = 64
    for n in range(1,n_max+1):
        #print("{0}\t{1}".format(n,compute_lower_bound(n)))
        print("{0}\t{1}".format(n, lower_bound(n)))
