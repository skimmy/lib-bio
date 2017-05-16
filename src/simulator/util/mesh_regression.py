import ub_regression as ubreg

import numpy as np

import sys
import math

def g(n, alpha, beta, gamma):
    return alpha*n + gamma*pow(n,beta)

def sum_of_squares_loss(tilde_e, model, weighted=False):
    squares_sum = 0
    weight = 1
    for p in tilde_e:
        residual = model(p[0]) - p[1]
        if (weighted):
            weight = p[3] / p[2]
        squares_sum += residual*residual * weight
    return squares_sum

def sum_of_modules_loss(tilde_e, model, wighted=False):
    module_sum = 0;
    weight = 1
    for p in tilde_e:
        if (weighted):
            weight = math.sqrt(p[3]) / math.sqrt(p[2])
        residual = model(p[0]) - p[1]
        module_sum += abs(residual) * weight
    return module_sum

def mesh_parameter_estimation(tilde_e, alpha_int, beta_int, gamma_int, loss_func, weighted=False):
    alpha_0 = alpha_int[0]
    alpha_1 = alpha_int[1]
    N_alpha = alpha_int[2]
    alpha_step = (alpha_1-alpha_0) / N_alpha

    beta_0 = beta_int[0]
    beta_1 = beta_int[1]
    N_beta = beta_int[2]
    beta_step = (beta_1-beta_0) / N_beta

    gamma_0 = gamma_int[0]
    gamma_1 = gamma_int[1]
    N_gamma = gamma_int[2]    
    gamma_step = (gamma_1-gamma_0) / N_gamma        

    mesh = np.zeros((N_alpha, N_beta, N_gamma))

    alpha = alpha_0
    for a in range(N_alpha):
        beta = beta_0
        for b in range(N_beta):
            gamma = gamma_0
            for c in range(N_gamma):
                mesh[a][b][c] = loss_func(tilde_e, lambda x: g(x,alpha,beta,gamma))
                gamma += gamma_step
                c +=1                
            beta += beta_step;
            b += 1            
        alpha += alpha_step
        a += 1
    return mesh

'''Returns a list with the coordinates of k minima from mesh'''
def get_minima_mesh_point(mesh, k=32):
    import heapq as hq
    h = []
    for a in range(mesh.shape[0]):
        for b in range(mesh.shape[1]):
            for c in range(mesh.shape[2]):
                
                hq.heappush(h, (mesh[a][b][c], (a, b, c)) )
    return hq.nsmallest(k,h)

    
def plot_mesh(mesh):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')


def plot_minima_curves(data, minima):

    ns = np.array([entry[0] for entry in data])
    tilde_e = np.array([entry[1] for entry in data])
    ns_g = ns[1:]
    tilde_g = np.array([2*tilde_e[i-1]-tilde_e[i] for i in range(1,len(ns))])
    print(tilde_g)
    
    import matplotlib.pyplot as plt
    fig, (ax, ax_g) = plt.subplots(1,2)

    # plot sample points
    ax.plot(ns, tilde_e, 'o')
    ax_g.plot(ns_g, tilde_g, 'o')

    avg_params = np.zeros(3)
    K = len(minima)    
    # plot all the |minima| curves
    for k in range(K):
        params = np.array(minima[k])
        avg_params = np.add(avg_params, params)
        y = np.array([g(n, params[0], params[1], params[2]) for n in ns])
        y_g = np.array([2*y[i-1] - y[i] for i in range(1,len(ns))])
        ax.plot(ns, y, color='r', alpha=0.25)
        ax_g.plot(ns_g, y_g, color='r', alpha=0.25)
    

    # plot the AVG(minima) curve
    avg_curve = [g(n,avg_params[0]/K, avg_params[1]/K, avg_params[2]/K) for n in ns]
    ax.plot(ns, avg_curve, 'r', alpha=0.75)

    plt.show()


if (__name__ == "__main__"):

    dp, tilde_en = ubreg.load_point(sys.argv[1])
    weighted = (sys.argv.count("-weight") > 0)
    mesh_edge = 16

    for rec_steps in range(1):

        alpha_i = (0.5, 0.52, mesh_edge)
        beta_i  = (0.34, 0.36, mesh_edge)
        gamma_i = (0, 1, mesh_edge)
        alpha_s = (alpha_i[1]-alpha_i[0]) / alpha_i[2]
        beta_s = (beta_i[1]-beta_i[0]) / beta_i[2]
        gamma_s = (gamma_i[1]-gamma_i[0]) / gamma_i[2]
    
        mesh = mesh_parameter_estimation(tilde_en, alpha_i, beta_i, gamma_i,
                                         sum_of_squares_loss, weighted)

        min_loss = np.amin(mesh)
        print(min_loss)

        minima = get_minima_mesh_point(mesh, 8)
        min_params = []
        for m in minima:            
            par = (alpha_i[0] + (m[1][0]*alpha_s), beta_i[0]+(m[1][1]*beta_s), gamma_i[0]+(m[1][2]*gamma_s))
            print("{0}\t{1}\t{2}".format(m[1], m[0], par))
            min_params.append(par)
        
    
    if (sys.argv.count("-plot") > 0):
        #plot_mesh(mesh)
        plot_minima_curves(tilde_en, min_params)
