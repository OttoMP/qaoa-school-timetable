# Import tools for running QAOA
from .qaoa import qaoa

import random

#import math tools
import numpy as np

# We import the tools to handle general Graphs
import networkx as nx
import matplotlib.pyplot as plt

# Import miscellaneous tools
from scipy.optimize import minimize
from math import floor, ceil
import os, psutil, datetime
import pandas as pd

def save_csv(data, nome_csv):
    data_points = pd.DataFrame(data, columns=['Expected Value', 'nfev', 'Beta0|Gamma|Beta'])
    data_points.to_csv(nome_csv, mode='a', header=False)
    return

def parameter_setting(gamma, beta, p):
    # -------------
    # Interpolation 
    # -------------
    next_gamma = [0]*(2*p)
    next_beta = [0]*(2*p)
    
    next_gamma[0] = gamma[0]
    next_beta[0] = beta[0]
    next_gamma[-1] = gamma[-1]
    next_beta[-1] = beta[-1]
    if p > 1:
        for i in range(1,2*p-1,2):
            next_gamma[i]   = (ceil(i/2)/p) * gamma[int(i/2)+1] + (floor(p-(i/2))/p) * gamma[int(i/2)]
            next_gamma[i+1] = (ceil(i/2)/p) * gamma[int(i/2)]   + (floor(p-(i/2))/p) * gamma[int(i/2)+1]
            next_beta[i]    = (ceil(i/2)/p) * beta[int(i/2)+1]  + (floor(p-(i/2))/p) * beta[int(i/2)]
            next_beta[i+1]  = (ceil(i/2)/p) * beta[int(i/2)]    + (floor(p-(i/2))/p) * beta[int(i/2)+1]
    
    return next_gamma, next_beta

def minimization_process_cobyla(goal_p, G, num_colors, school, cost_function):
    iterations = 10 # Number of independent runs
    
    local_optima_param = []
    # --------------------------
    # COBYLA Optimization
    # --------------------------
    for i in range(iterations):
        p = 1          # Start value of p
        while p <= goal_p:
            with open(f'timing/{school}_{p}_{i}.txt', 'w') as file:
                pass
            epsilon = 1e-6
            maxfev = 100
            if p == 8:
                epsilon = 1e-4
                maxfev = 50
            qaoa_args = p, G, num_colors, epsilon, cost_function, school, i
            print("Running minimization process with p-value", p)
            # --------------------------
            # Initializing QAOA Parameters 
            # --------------------------
            if p > 1:
                # Extracting previous local optima
                beta0 = local_optima_param[0]
                new_local_optima_param = np.delete(local_optima_param, 0)
                middle = int(len(local_optima_param)/2)
                p_gamma = new_local_optima_param[:middle] # Previous gamma
                p_beta = new_local_optima_param[middle:]  # Previous beta
                
                # Parameter setting strategy
                gamma, beta = parameter_setting(p_gamma, p_beta, int(p/2))
            else:
                beta0 = random.uniform(0, np.pi)
                gamma = [random.uniform(0, 2*np.pi)]
                beta  = [random.uniform(0, np.pi)]
            #print("Using Following parameters:")
            #print("Beta0:", beta0)
            #print("Gamma:", gamma)
            #print("Beta:", beta)
            qaoa_par = [beta0]+gamma+beta
            
            # Construct parameters bounds in the form of constraints
            beta0_bounds = [[0, np.pi]]
            beta_bounds = [[0, np.pi]]*p
            gamma_bounds = [[0, 2*np.pi]]*p
            bounds = beta0_bounds+gamma_bounds+beta_bounds
            cons = []
            for factor in range(len(bounds)):
                lower, upper = bounds[factor]
                l = {'type': 'ineq',
                    'fun': lambda x, lb=lower, i=factor: x[i] - lb}
                u = {'type': 'ineq',
                    'fun': lambda x, ub=upper, i=factor: ub - x[i]}
                cons.append(l)
                cons.append(u)
            
            #print("\nMemory Usage", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
            print("Minimizing function using COBYLA")
            print("Current Time:-", datetime.datetime.now())
            res = minimize(qaoa, qaoa_par, args=qaoa_args, method='COBYLA',
                    constraints=cons, options={'disp': False, 'maxiter': maxfev})
            print(res)
            print("Current Time:-", datetime.datetime.now())
            #print("Memory Usage", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
            print("Saving Results\n")
            save_csv([[res['fun'], res['nfev'], res['x']]], f"results/{school}_{p}_{i}.csv" )
            local_optima_param = res['x']
            
            # Preparing next p-value
            p = p*2
