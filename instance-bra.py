# Import tools for running QAOA
import random
import cma
from ket import *

#import math tools
import numpy as np

# We import the tools to handle general Graphs
import networkx as nx

# Import miscellaneous tools
from scipy.optimize import minimize
from math import floor, ceil
import sys, os, psutil, datetime
from xml_parser import parseXML
from itertools import combinations
import pprint as pp
import pandas as pd

def save_csv(data, nome_csv):
    data_points = pd.DataFrame(data, columns=['Expected Value', 'Beta0|Gamma|Beta'])
    data_points.to_csv(nome_csv, mode='a', header=False)

    return

def cost_function_bra(G):
    C = 0

    # Find times for each teacher
    times_T1 = []
    times_T2 = []
    times_T3 = []
    times_T4 = []
    times_T5 = []
    times_T6 = []
    times_T7 = []
    times_T8 = []

    for node in G.nodes:
        if node.find("T1") != -1:
            times_T1.append(G.nodes[node]['color'])
        if node.find("T2") != -1:
            times_T2.append(G.nodes[node]['color'])
        if node.find("T3") != -1:
            times_T3.append(G.nodes[node]['color'])
        if node.find("T4") != -1:
            times_T4.append(G.nodes[node]['color'])
        if node.find("T5") != -1:
            times_T5.append(G.nodes[node]['color'])
        if node.find("T6") != -1:
            times_T6.append(G.nodes[node]['color'])
        if node.find("T7") != -1:
            times_T7.append(G.nodes[node]['color'])
        if node.find("T8") != -1:
            times_T8.append(G.nodes[node]['color'])

    # No IDLE time for teachers constraint Soft 3
    # Counting IDLE time
    for i in times_T1[:-1]:
        for ii in times_T1[i:]:
            # If in the same day and is not consecutive
            if int(i/3) == int(ii/3) and abs(i-ii) == 2:
                C += 3
        
    # Not more than 2 days with lessons for teacher Soft 9
    # Not more than 3 days with lessons for teacher Soft 9 - only T6
    # Count number of days from teacher
    days_T1 = len(set([int(time/3) for time in times_T1]))
    days_T2 = len(set([int(time/3) for time in times_T2]))
    days_T3 = len(set([int(time/3) for time in times_T3]))
    days_T4 = len(set([int(time/3) for time in times_T4]))
    days_T5 = len(set([int(time/3) for time in times_T5]))
    days_T6 = len(set([int(time/3) for time in times_T6]))
    days_T7 = len(set([int(time/3) for time in times_T7]))
    days_T8 = len(set([int(time/3) for time in times_T8]))

    days_list = [days_T1, days_T2, days_T3, days_T4, days_T5, days_T7, days_T8]
    for days in days_list:
        if days > 2:
            C += 9
    if days_T6 > 3:
        C += 9

    return C

def partial_mixer(qc, neighbour, ancilla, target, beta):
    def outer():
        if neighbour == None:
            X(ancilla)
        else:
            with around(X, neighbour):
                ctrl(neighbour, X, ancilla)

    with around(outer):
        with around([H, ctrl(0, X, target=1)], target):
            ctrl(ancilla, RZ, 2*beta, target[1])
        
        with around([RX(-np.pi/2), ctrl(0, X, target=1)], target):
            ctrl(ancilla, RZ, 2*beta, target[1])

def neighbourhood(G, num_colors, node, color, list_nodes):
    neighbours = list(G[node])
    neighbours_index = [list_nodes.index(neigh) for neigh in neighbours]

    neighbours_color_qubit = [color+(num_colors*u) for u in neighbours_index]

    return neighbours_color_qubit

# Apply the partial mixer for each pair of colors of each node
def mixer(qc, G, beta, num_nodes, num_colors):
    list_nodes = list(G.nodes())
    for u, node in enumerate(G.nodes):
        for i in range(num_colors):
            neighbours_i = neighbourhood(G, num_colors, node, i, list_nodes)
            for j in range(num_colors):
                if i < j:
                    neighbours_j = neighbourhood(G, num_colors, node, j, list_nodes)
                    neighbours = neighbours_i+neighbours_j

                    if neighbours == []:
                        q_neighbours = None
                    else:
                        q_neighbours = qc[neighbours[0]]
                        for neigh in neighbours[1:]:
                            q_neighbours = q_neighbours | qc[neigh]
                    partial_mixer(
                            qc,
                            q_neighbours,
                            qc[num_nodes*num_colors+u],
                            qc[i+(num_colors*u)]|qc[j+(num_colors*u)],
                            beta)

def phase_separator(qc, gamma, num_nodes, num_colors):
    #qubits2color(qc, num_nodes, num_colors)
    for node in range(num_colors*num_nodes):
        X(qc[node])
    for k in range(num_colors):
        qubits = [node*num_colors+k for node in range(num_nodes)]
        control = qc[qubits[0]]
        for qub in qubits[1:-1]:
            control = control | qc[qub]
        target = qc[qubits[-1]]
        ctrl(control, RZ, 2*gamma, target)
    for node in range(num_colors*num_nodes):
        X(qc[node])
    #color2qubits(qc, num_nodes, num_colors)

def qaoa_min_graph_coloring(p, G, num_nodes, num_colors, beta0, gamma, beta):
    # --------------------------
    # Initializing qubits
    # --------------------------
    qc = quant((num_nodes*num_colors) + num_nodes)
    # --------------------------
    # Initial state preparation
    # --------------------------
    coloring = [G.nodes[node]['color'] for node in G.nodes]
    for i, color in enumerate(coloring):
        X(qc[(i*num_colors)+color])

    # --------------------------
    # Alternate application of operators
    # --------------------------
    mixer(qc, G, beta0, num_nodes, num_colors) # Mixer 0
    for step in range(p):
        phase_separator(qc, gamma[step], num_nodes, num_colors)
        mixer(qc, G, beta[step], num_nodes, num_colors)

    # --------------------------
    # Measurement
    # --------------------------
    #result = measure(qc).get()
    return dump(qc)

def qaoa(par, p, initial_G, num_colors):
    # --------------------------
    # Unpacking QAOA parameters
    # --------------------------
    beta0 = par[0]
    new_par = np.delete(par, 0)
    middle = int(len(par)/2)
    gamma = new_par[:middle]
    beta = new_par[middle:]
    num_nodes = initial_G.number_of_nodes()

    # --------------------------
    # Verifying Parameters
    # --------------------------
    #print("Using Following parameters:\nBeta0:", beta0, "\nGamma:", gamma, "\nBeta:", beta)

    # --------------------------
    # Running QAOA on simulator
    # --------------------------
    G = nx.Graph()
    G.add_nodes_from(initial_G)
    G.add_edges_from(initial_G.edges)
    initial_coloring = [initial_G.nodes[node]['color'] for node in initial_G.nodes]
    color_graph_from_coloring(G, initial_coloring)
    
    result = qaoa_min_graph_coloring(p, initial_G, num_nodes, num_colors, beta0, gamma, beta)
    #print("Number of States", len(result.get_states()))
    #print("State Vector", result.show('b6:b6:b6:b6:b6:b6'))

    # --------------------------
    # Counting resulting states
    # --------------------------
    counts = {} # Dictionary for keeping the results of the simulation
    for i in result.get_states():
        binary = np.binary_repr(i, width=(num_nodes*num_colors)+num_nodes)
        counts[binary] = int(2**20*result.probability(i))

    # --------------------------
    # Evaluate the data from the simulator
    # --------------------------
    avr_function_value = 0
    min_function_value = [0, np.inf]

    for sample in list(counts.keys()):
        if counts[sample] > 0:
            # use sampled bit string x to compute f(x)
            x       = [int(num) for num in list(sample)]
            
            # Coloring Graph with counts[sample]
            coloring = []
            for i in range(len(G)):
                for pos, char in enumerate(x[i*num_colors:(i*num_colors+num_colors)]):
                    if int(char):
                        coloring.append(pos)
            color_graph_from_coloring(G, coloring)

            fx = cost_function_bra(G)

            # compute the expectation value and energy distribution
            avr_function_value = avr_function_value + counts[sample]*fx

            # save best bit string
            if( min_function_value[1] > fx):
                min_function_value[0] = sample
                min_function_value[1] = fx

    expectation_value = avr_function_value/sum(counts.values())

    return expectation_value

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

def minimization_process_cobyla(goal_p, G, num_colors, school):
    iterations = 50 # Number of independent runs
    
    local_optima_param = []
    # --------------------------
    # COBYLA Optimization
    # --------------------------
    for i in range(iterations):
        p = 1          # Start value of p
        qaoa_args = p, G, num_colors
        while p <= goal_p:
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
                    constraints=cons, options={'disp': False})
            print(res)
            print("Current Time:-", datetime.datetime.now())
            #print("Memory Usage", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
            print("Saving Results\n")
            save_csv([[res['fun'], res['x']]], "results/cobyla/"+school+"p"+str(p)+".csv" )
            local_optima_param = res['x']
            
            # Preparing next p-value
            p = p*2

def minimization_process_cma(goal_p, G, num_colors, school): 
    # --------------------------
    # CMA-ES Optimization 
    # --------------------------
    local_optima_param = []
    p = 1          # Start value of p
    while p <= goal_p:
        print("Running minimization process with p-value", p)
        #print("\nMemory Usage", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
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

        # Settings parameters bounds
        lower_bounds = [0] * ((2*p)+1)
        upper_bounds_beta0 = [np.pi]
        upper_bounds_gamma = [2*np.pi]*p
        upper_bounds_beta  = [np.pi]*p
        upper_bounds = upper_bounds_beta0+upper_bounds_gamma+upper_bounds_beta 
        opts = {'bounds' : [lower_bounds, upper_bounds], 'maxiter': 300, } #'maxfevals': 300}
        sigma0 = 0.3*(2*np.pi)
        print("Initial Step =", sigma0)
        
        es = cma.CMAEvolutionStrategy(qaoa_par, sigma0, opts)
        while not es.stop():
            solutions = es.ask()
            function_values = [qaoa(s, p, G, num_colors) for s in solutions]
            es.tell(solutions, function_values)
            res = es.result
            #print("Saving Results")
            #save_csv([[res[1], res[0]]], "results/cma/"+school+"p"+str(p)+".csv" )
            es.disp()
        print("---------------------------")
        es.result_pretty()
        res = es.result
        print("---------------------------")
        print("Optimal Result", res[0])
        print("Respective Function Value", res[1])
        print("Respective Function Evaluations", res[2])
        print("Overall Function Evaluations", res[3])
        print("Overall Iterations", res[4])
        print("Mean Result", res[5])
        print("Standard Deviation Final Sample", res[6])
        
        #print("Memory Usage", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
        print("Saving Final Results")
        print("---------------------------\n")
        save_csv([[res[1], res[0]]], "results/cma/"+school+"p"+str(p)+".csv" )
        local_optima_param = res[0]

        # Preparing next p-value
        p = p*2

def color_graph_from_coloring(graph, coloring):
    for index,node in enumerate(graph.nodes):
        graph.nodes[node]['color'] = coloring[index]

    return

def color_graph_from_num(graph, num_color):
    color_index = 0
    node_list = list(graph.nodes)
    not_allowed_color = []

    # Mark all the vertices as not visited
    visited = {x: False for x in graph.nodes}
    for node in graph.nodes:
        if graph.nodes[node]['color'] != None:
            visited[node] = True

    # Create a queue for BFS
    queue = []

    # Mark the source node as
    # visited and enqueue it
    for u in node_list:
        success = True
        queue.append(u)
        visited[u] = True

        while queue:
            # Dequeue a vertex from queue and color it
            source = queue.pop(0)
            not_allowed_color = {graph.nodes[neighbour]['color'] for neighbour in graph[source]
                                    if (graph.nodes[neighbour]['color'] != None) }
            if len(not_allowed_color) == num_color:
                success = False
                visited = {x: False for x in graph.nodes}
                break

            while color_index in not_allowed_color:
                color_index = (color_index+1)%num_color
            graph.nodes[source]['color'] = color_index
            not_allowed_color = set()

            # Get all adjacent vertices of the
            # dequeued vertex s. If a adjacent
            # has not been visited, then mark it
            # visited and enqueue it
            for i in graph[source]:
                if visited[i] == False:
                    queue.append(i)
                    visited[i] = True
        if success:
            break

    return

def create_graph_from_events(events):
    G = nx.Graph()
    
    week = [
            "Mo_1",
            "Mo_2",
            "Mo_3",
            "Tu_1",
            "Tu_2",
            "Tu_3",
            "We_1",
            "We_2",
            "We_3",
            "Th_1",
            "Th_2",
            "Th_3",
            "Fr_1",
            "Fr_2",
            "Fr_3",
    ]
    comb = combinations(week, 2)
    G.add_nodes_from(week)
    G.add_nodes_from([(day, {'color' : c}) for c, day in enumerate(week)])
    for i in comb:
        G.add_edge(i[0], i[1])
    
    # Adding events to graph
    # No clash constraints Hard 1
    G.add_nodes_from([(event['Id'], {'color' : None}) for event in events])
    comb = combinations(events, 2)
    for i in comb:
        res0 = set(i[0]['Resources'])
        res1 = i[1]['Resources']
        intersection = [value for value in res0 if value in res1]
        if intersection:
            G.add_edge(i[0]['Id'], i[1]['Id'])
            
    # Forbidden time constraint Hard 1
    #  T1 We
    G.add_edge("T1-S1", "We_1")
    G.add_edge("T1-S1", "We_2")
    G.add_edge("T1-S1", "We_3")
    G.add_edge("T1-S2", "We_1")
    G.add_edge("T1-S2", "We_2")
    G.add_edge("T1-S2", "We_3")
    G.add_edge("T1-S3", "We_1")
    G.add_edge("T1-S3", "We_2")
    G.add_edge("T1-S3", "We_3")
    #  T2 Fr
    G.add_edge("T2-S1", "Fr_1")
    G.add_edge("T2-S1", "Fr_2")
    G.add_edge("T2-S1", "Fr_3")
    G.add_edge("T2-S2", "Fr_1")
    G.add_edge("T2-S2", "Fr_2")
    G.add_edge("T2-S2", "Fr_3")
    #  T3 Th
    G.add_edge("T3-S1", "Th_1")
    G.add_edge("T3-S1", "Th_2")
    G.add_edge("T3-S1", "Th_3")
    G.add_edge("T3-S2", "Th_1")
    G.add_edge("T3-S2", "Th_2")
    G.add_edge("T3-S2", "Th_3")
    G.add_edge("T3-S3", "Th_1")
    G.add_edge("T3-S3", "Th_2")
    G.add_edge("T3-S3", "Th_3")
    #  T4 Tu
    G.add_edge("T4-S1", "Tu_1")
    G.add_edge("T4-S1", "Tu_2")
    G.add_edge("T4-S1", "Tu_3")
    G.add_edge("T4-S2", "Tu_1")
    G.add_edge("T4-S2", "Tu_2")
    G.add_edge("T4-S2", "Tu_3")
    G.add_edge("T4-S3", "Tu_1")
    G.add_edge("T4-S3", "Tu_2")
    G.add_edge("T4-S3", "Tu_3")
    #  T5 Mo
    G.add_edge("T5-S2", "Mo_1")
    G.add_edge("T5-S2", "Mo_2")
    G.add_edge("T5-S2", "Mo_3")
    G.add_edge("T5-S3", "Mo_1")
    G.add_edge("T5-S3", "Mo_2")
    G.add_edge("T5-S3", "Mo_3")
    #  T6 We
    G.add_edge("T6-S1", "We_1")
    G.add_edge("T6-S1", "We_2")
    G.add_edge("T6-S1", "We_3")
    G.add_edge("T6-S2", "We_1")
    G.add_edge("T6-S2", "We_2")
    G.add_edge("T6-S2", "We_3")
    G.add_edge("T6-S3", "We_1")
    G.add_edge("T6-S3", "We_2")
    G.add_edge("T6-S3", "We_3")
    #  T7 Tu
    G.add_edge("T7-S1", "Tu_1")
    G.add_edge("T7-S1", "Tu_2")
    G.add_edge("T7-S1", "Tu_3")
    G.add_edge("T7-S3", "Tu_1")
    G.add_edge("T7-S3", "Tu_2")
    G.add_edge("T7-S3", "Tu_3")
    #  T8 Th
    G.add_edge("T8-S1", "Th_1")
    G.add_edge("T8-S1", "Th_2")
    G.add_edge("T8-S1", "Th_3")
    G.add_edge("T8-S2", "Th_1")
    G.add_edge("T8-S2", "Th_2")
    G.add_edge("T8-S2", "Th_3")
    G.add_edge("T8-S3", "Th_1")
    G.add_edge("T8-S3", "Th_2")
    G.add_edge("T8-S3", "Th_3")

    return G

def main():
    print("Starting program\n")

    # --------------------------
    # School Instances
    # --------------------------
    school = "Bra"

    # --------------------------
    # Parse XML file
    # --------------------------
    events = parseXML('dataset/bra-instance01.xml')

    # --------------------------
    #  Preparing Conflict Graph
    # --------------------------
    G = create_graph_from_events(events)
    
    print("--------------------------")
    print("Graph information\n")
    
    print("Nodes = ", G.nodes)
    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nPre-coloring", coloring)

    degree = [deg for (node, deg) in G.degree()]
    print("\nDegree of each node", degree)

    # --------------------------
    #  Coloring Conflict Graph
    # --------------------------
    num_colors = 15 
    color_graph_from_num(G, num_colors)
    
    coloring = [G.nodes[node]['color'] for node in G.nodes]
    
    print("\nNumber of colors", num_colors)
    print("\nInitial coloring", coloring)

    initial_function_value = cost_function_bra(G)
    print("\nInitial Function Value", initial_function_value)

    # ---------------------------
    # Verifying Graph consistency
    #----------------------------
    print("----------------------------")
    print("Verifying Graph consistency")
    for i in G.nodes:
        print("\nNode",i,"Color", G.nodes[i]['color'])
        color_and_neighbour = [(neighbour, G.nodes[neighbour]['color']) for neighbour in G[i]]
        print("Neighbours | Color")
        for pair in color_and_neighbour:
            print(pair)
            
    #----------------------------
    # Starting QAOA
    #----------------------------
    print("----------------------------")
    print("Running QAOA")
    num_nodes = G.number_of_nodes()
    number_of_qubits = num_nodes*num_colors+num_nodes
    print("Necessary number of qubits: ", number_of_qubits)

    # QAOA parameter
    goal_p = 8

    # Minimizing Example Bra
    minimization_process_cobyla(goal_p, G, num_colors, school)
    minimization_process_cma(goal_p, G, num_colors, school)

    print("Program End")
    print("----------------------------")

if __name__ == '__main__':
    main()