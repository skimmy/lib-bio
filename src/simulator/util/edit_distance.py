def edit_distance(s1, s2, verbose=False):
    import numpy as np
    n1 = len(s1)+1
    n2 = len(s2)+1
    dp_matrix = np.zeros([n1, n2])
    for i in range(n1):
        dp_matrix[i,0] = i
    for j in range(n2):
        dp_matrix[0,j] = j

    for i in range(1,n1):
        for j in range(1,n2):
            delta = 1
            if (s1[i-1] == s2[j-1]):
                delta = 0
            dp_matrix[i,j] = min(dp_matrix[i-1,j-1]+delta,
                                 dp_matrix[i-1,j]+1,
                                 dp_matrix[i,j-1]+1)        
    if (verbose):
        print dp_matrix    
    return int(dp_matrix[n1-1,n2-1])

if (__name__ == "__main__"):
    import sys
    s1 = sys.argv[1]
    s2 = sys.argv[2]
    verbose=False
    ED = edit_distance(s1, s2, verbose)
    print("{0}\t{1}\t{2}".format(s1, s2, ED))
