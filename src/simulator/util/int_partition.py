import numpy as np
from scipy.special import binom

def rec_part(n,k):
    if (n==0):
        return []
    if (k==1):
        return [ [n] ]
    P = []
    for p in range(1,n+1):
        rs  = rec_part(n-p, k-1)
        for r in rs:
            P.append([p] + r)            
    return P

def count_positive(v):
    x = 0
    for i in v:
        if (i>0):
            x +=1
    return x

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    sigma_max = 4
    Sigma = [str(x) for x in range(sigma_max+1)]
    SigmaR = np.arange(sigma_max,0,-1)
    BaseChoices = np.arange(1,sigma_max+1)
    cumsum = 0
    totdistinct = 0
    out_file = open("./int_part_" + str(n) + ".txt", "w")
    out_file_str = open("./part_strings_" + str(n) + ".txt", "w")
    for sigma in range(1,sigma_max+1):        
        parts = rec_part(n,sigma)        
        for P in parts:
            totdistinct += 1
            P.extend([0]*(sigma_max-len(P)))
            out_file.write(",".join(str(x) for x in P) + "\n")
            S = ""
            for i in range(sigma_max):
                S += Sigma[i]*P[i]
            B = np.prod(SigmaR[:sigma])
            M = np.power(BaseChoices, (np.array(P)-1).clip(min=0)).prod()
            out_file_str.write(S + "\n")
            cumsum += B * M
            print("{0}\t{1}x{2}\t{3}".format(P, B, M, S))
    print(cumsum)
    print(np.power(sigma_max, n))
    out_file.close()
    out_file_str.close()
