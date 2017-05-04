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

'''Returns the smallest integer greater than or equal to x and divisible by n'''
def get_ceil_divisible_by(x, n):
    x_ceil = int(math.ceil(x))
    k = 1
    while(k < x_ceil):
        k += n
    return k

def plot_curves(params_dict, save_file=False):
    import matplotlib.pyplot as plt

    # dictionary items
    # (beta, gamma, err, datapoints, logpoints)
    plt_colors = ['b', 'r', 'g', 'k', 'y']
    color_i = 0
    fig, ax = plt.subplots()
    #ax.set_autoscale_on(False)
    max_val = 0
    max_n = 0
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
            logn = math.log(n,2)
            ns.append(n)
            exp_points.append(dpoints[i][1])
            y.append(pow(2, gamma+beta*logn))            
        style_curve = plt_colors[color_i] + '-'
        style_point = plt_colors[color_i] + 'o'
        ax.plot(ns, exp_points, style_point)
        ax.plot(ns, y, style_curve)
        color_i += 1
        max_n = max(max_n,max(ns))
        max_val = max(max_val, max(max(y), max(exp_points) ))
        
        
    plt.title("Sample points and regression curve")
    plt.xlabel("n")
    plt.ylabel("2g(n/2) - g(n)")
    plt.xlim(xmin=0)
    print(max_val)
    y_max = get_ceil_divisible_by(max_val,5)
    plt.ylim([0, y_max])
    plt.grid(True)
    if (save_file):
        plt.savefig("plot.pdf")
    else:
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

def test(log_p, g, b):
    pass

if (__name__ == "__main__"):
    print("")
    file_names = sys.argv[1].split(',')
    params_dict = dict()
    # for each given file
    for file_name in file_names:
        # load and transform into log space
        datapoints, en = load_point(file_name)
        logpoints = loglog_scale(datapoints)
        # regression
        beta_p, gamma_p, err = regression(logpoints)
        params_dict[file_name] = (beta_p, gamma_p, err, datapoints, logpoints)
        print("File: {0}\tm = {1}".format(file_name, len(en)))
        print("\n".join([str(x[1]) for x in datapoints]))
        print("")
        print("Gamma'  = {0}\nBeta'   = {1}\nError  = {2}".format(gamma_p, beta_p, err))        
        # if test is mode is on
        if (sys.argv.count("-test") > 0):
            test(logpoints, gamma_p, beta_p)
            exit(0)
    # plots all the curves and sample points
    if (sys.argv.count("-plot") > 0):
        save= False
        if(sys.argv.count("-save") > 0):
            save = True
        plot_curves(params_dict, save)
