import sys
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
        d_bar = min(q-1,D)
        for d in range(d_bar):
            count += (spec.binom(q-1,d) * spec.binom(D-1,d)
                      * spec.binom(D-1, q-d-1))
        count *= pow(Sigma-1, D)
        
        # ending ins
        i_bar = min(q-1, D)
        for i in range(i_bar):
            part_count = 0
            for j in range(1,D-i):
                part_count += (spec.binom(D-j-1,i) * spec.binom(D-1,q-i-1) *
                               pow(Sigma,j) * pow(Sigma-1,D-j))
            count += part_count * spec.binom(q-1,i)
            
    return count

##############################################################################
#                                                                            #
#                          NEW BOUND FUNCTIONS                               #
#                                                                            #
##############################################################################

def count_canonical_annotated_path_r_D_q(n, r, D, q, Sigma):
    count = 0
    # sum of fd*fn
    for delta_s in range(2):
        for delta_e in range(2):
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


def bound_for_hulls(n, Sigma):
    up_r = []
    for r in range(1,n+1):
        bound_r = count_canonical_annotated_path_r(n, r, Sigma)
        up_r.append((r,bound_r))
    return up_r
    


def loose_bound(n, Sigma=4):
    hulls_bound = list([(0,1)]) + bound_for_hulls(n,Sigma)
    lb = sum([entry[0]*entry[1] for entry in hulls_bound])
    volumes_bound = [0]*len(hulls_bound)
    volumes_bound[0] = hulls_bound[0][1]
    for i in range(1,len(hulls_bound)):
        volumes_bound[i] = volumes_bound[i-1] + hulls_bound[i][1]
    ub = n*pow(Sigma,n) - sum([volumes_bound[r] for r in range(1,n)])
    return (lb , ub, hulls_bound)

def bound(n, Sigma=4):
    lb = 0
    ub = 0
    hulls_bound = [(i,0) for i in range(n+1)]

    
    remaining_strings = pow(Sigma,n) - 1
    r = 1    
    while (remaining_strings > 0 and r <= n):
        hull_count = count_canonical_annotated_path_r(n,r, Sigma)     
        pay = min(hull_count, remaining_strings)
        lb += pay*r
        remaining_strings = remaining_strings - hull_count
        hulls_bound[r] = (r,lb)
        r += 1
    
    return (lb, ub, hulls_bound)

##############################################################################
#                                                                            #
#                                 MAIN                                       #
#                                                                            #
##############################################################################

if __name__ == "__main__":
    verbosity = 0
    if (len(sys.argv) > 1):
        verbosity = int(sys.argv[1])
    
    n_max = 128
    

    for n in range(1,n_max+1):
        norm_const = float(n * pow(4,n))
        (lb, ub, hulls_bound)  = bound(n)
        print("      n = {0}      ".format(n))
        if (verbosity > 0):
            print("  r\tS_r\n----------------")
            for entry in hulls_bound:
                print("  {0}\t{1}".format(entry[0],entry[1]))
        print("\nLB\t{1}\n".format(n, lb / norm_const))
