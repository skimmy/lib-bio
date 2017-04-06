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
import math
from scipy import stats

def g(n, gamma, beta):
    return gamma*pow(n, beta)

def plot_curves(params_dict):
    import matplotlib.pyplot as plt

    # dictionary items
    # (beta, gamma, err, datapoints, logpoints)
    plt_colors = ['b', 'r', 'g', 'k', 'y']
    color_i = 0
    fig, ax = plt.subplots()
    for k,v in params_dict.items():
        name = k
        beta = v[0]
        gamma = v[1]
        dpoints = v[3]
        ns = []
        exp_points = []
        y = []
        for i in range(1,len(dpoints)):
            n = dpoints[i][0]
            ns.append(n)
            exp_points.append(dpoints[i][1])
            y.append(2*g(n/2.0, gamma, beta) - g(n, gamma, beta))

        style_curve = plt_colors[color_i] + '-'
        style_point = plt_colors[color_i] + 'o'
        ax.plot(ns, exp_points, style_point)
        ax.plot(ns, y, style_curve)
        color_i += 1
                
    plt.show()


def load_point(file_name):
    in_file = open(file_name)
    # load e(ni) for i=0,..,m
    e_n = []
    for line in in_file:
        if (line[0] == '#'):
            continue
        (str_x, str_y) = line.strip().split()
        e_n.append((float(str_x), float(str_y)))        
    in_file.close()
    # transform
    datapoints = []
    for i in range(1,len(e_n)):
        n = e_n[i][0]
        en = e_n[i][1]
        en_half = e_n[i-1][1]
        datapoints.append( (float(n), 2.0*en_half - en) )
    return datapoints, e_n

def loglog_scale(line_points, base=2):
    log_points = [( math.log(x[0],base), math.log(x[1],base) ) for x in line_points]    
    return log_points

def regression(points):
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    slope, intercept, r_v, p_v, std_err = stats.linregress(x,y)
    return slope, intercept, std_err

if (__name__ == "__main__"):
    print()
    file_names = sys.argv[1].split(',')
    params_dict = dict()
    for file_name in file_names:
        datapoints, en = load_point(file_name)
        print (datapoints)
        logpoints = loglog_scale(datapoints)
        beta_p, gamma_p, err = regression(logpoints)
        beta = beta_p
        gamma =  pow(2,gamma_p) / (pow(2,beta+1) - 1)
        params_dict[file_name] = (beta, gamma, err, datapoints, logpoints)
        print("File: {0}\tm = {1}".format(file_name, len(en)))
        print()
        print("Gamma'  = {0}\nBeta'   = {1}\nError  = {2}".format(gamma_p, beta_p, err))
        print()
        print("Gamma  = {0}\nBeta   = {1}".format(gamma, beta))
        print()
    if (sys.argv.count("-plot") > 0):
        plot_curves(params_dict)
