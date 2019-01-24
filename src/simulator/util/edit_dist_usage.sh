# This file contains examples of usage of 'simulator' for tasks
# related to edit distance.

# !!! WORK IN PROGRESS !!!

################## EXHASUTIVE ####################

# Computes exact expected edit distance

# FLAGS

# -O 6 -->  Edit Distance Task
# -f 1 -->  Exhaustive calculation

./simulator.out -O 6 -f 1 -N 4 -t 2


# Computes distribution of eccentricity for all strings

# FLAGS

# -O 6 -->  Edit Distance Task
# -B 4 -->  Subtask
# -f 1 -->  Flags
# -t t -->  Number of threads used


./simulator.out -O 6 -B 4 -f 1 -N 4

# PARAMETERS

# -N   length of the strings


#################### SAMPLING ####################


#################### ESTIMATE ####################
# Estimate of ~g(n) through bandiwdth algorithm


# FLAGS

# -O 6 -->  Edit Distance Task
# -f 64 --> EDIT_DISTANCE_DIFF_BOUND_ERROR
#           Simulation stops when |~g - ~err| < f(epsilon, delta)
#           f(e,d) is a function s.t. P[err>e] < delta
#           Simulation may stop if k_max iterations have been performed
# -A x  --> Approximation type
#           -A 1   Estimate optimal bandwith once for all
#           -A 2   At each step estimate minimum bandwidth without
#                  variation w.r.t. the previous one.


# PARAMTERS

# -N    length of strings
# -k    max iteration
# -P    epsilon
# -c    delta


# EXAMPLE

./simulator.out -O 6 -N 8192 -A 2 -f 64 -P 0.05 -c 2 -k 50000


# NOTE

# The function called is
#    edit::difference_estimate_adaptive
# from
#    include/edit_estimates.hpp
