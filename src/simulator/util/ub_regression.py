# Performs regression for constants beta and gamma representing the
# parameters of the 'second order' function of alpha.
# Regression is perfomed on the input file which must have the format
#
#    n0   e(n0)
#    n1   e(n1)
#      ...
#
# scipy library is used to perform actual regression see
#   scipy.stats.lingress

import sys

if (__name__ == "__main__"):
    file_name = sys.argv[1]
    in_file = open(file_name)
    datapoints = []
    for line in in_file:
        (str_x, str_y) = line.strip().split()
        datapoints.append((float(str_x), float(str_y)))        
    in_file.close()
    print(datapoints)
