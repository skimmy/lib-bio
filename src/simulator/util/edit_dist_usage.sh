# This file contains examples of usage of 'simulator' for tasks
# related to edit distance.

# !!! WORK IN PROGRESS !!!

#################### SAMPLING ####################


#################### ESTIMATE ####################
# Estimate of ~g(n) through approximated algorithm


# FLAGS

# -O 6 -->  Edit Distance Task
# -f 64 --> EDIT_DISTANCE_DIFF_BOUND_ERROR
#           Simulation stops when |~g - ~err| < f(epsilon, delta)
#           f(e,d) is a function s.t. P[err>e] < delta
#           Simulation may stop if k_max iterations have been performed
# -A x  --> Use approximated edit distance algorithm when x > 0


# PARAMTERS

# -N    length of strings
# -k    max iteration
# -P    epsilon
# -c    delta


# EXAMPLE

./simulator.out -O 6 -N 8192 -A 1 -f 64 -P 0.05 -c 2 -k 50000


# NOTE

# The function called is
#    edit::difference_estimate_adaptive
# from
#    include/edit_estimates.hpp
