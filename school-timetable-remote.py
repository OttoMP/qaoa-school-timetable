# Import tools for running QAOA
import random
import cma
from scipy.optimize import minimize
from ket import *
from ket.lib import swap

#import math tools
import numpy as np

# We import the tools to handle general Graphs
import networkx as nx

# Import miscellaneous tools
from xml_parser import parseXML
from itertools import combinations
import pprint as pp
from binarytree import build

# We import plotting tools
import pandas as pd
import matplotlib.pyplot as plt

# Compute the value of the cost function
def cost_function_timetable(x, G, num_colors, list_students):

    coloring = []
    for i in range(len(G)):
        for pos, char in enumerate(x[i*num_colors:(i*num_colors+num_colors)]):
            if int(char):
                coloring.append(pos)

    color_graph_coloring(G, coloring)

    C = 0
    lectures = G.nodes
    for student in list_students:
        lecture_student = [(j,k,l) for (j,k,l) in lectures if l == student]
        timeslots = [G.nodes[key]['color'] for key in list(lecture_student)]
        timeslots.sort()
        new_timeslots = [x%num_colors for x in timeslots]

        # Students shouldn't have lectures at the last time of the day
        for time in new_timeslots:
            if time % num_colors == num_colors-1:
                C += 1

        day = []
        last_lecture = -np.inf
        for time in new_timeslots:
            if last_lecture >= time:
                # Students shouldn't have only one lecture in a day
                if len(day) == 1:
                    C += 1
                #Students shouldn't have many consecutive lectures
                else:
                    for index, lect in enumerate(day):
                        if index+1 < len(day) and day[index+1] == lect+1:
                            C += 1
                day = []

            last_lecture = time
            day.append(last_lecture)
        for index, lect in enumerate(day):
            if index+1 < len(day) and day[index+1] == lect+1:
                C += 1

    return C

def initial_cost_function_den_25pts(G):
    C = 0
    
    if G.nodes["Event1"]['color'] != 0:
        C += 1
    if G.nodes["Event2"]['color'] != 2:
        C += 1
    if G.nodes["Event3"]['color'] != 3:
        C += 1
    if G.nodes["Event4"]['color'] != 1:
        C += 1
    if G.nodes["Event5"]['color'] != 2:
        C += 1
    if G.nodes["Event6"]['color'] != 3:
        C += 1
    if G.nodes["Event7"]['color'] != 3:
        C += 1
    if G.nodes["Event8"]['color'] != 2:
        C += 1
    if G.nodes["Event9"]['color'] != 0:
        C += 1
    if G.nodes["Event10"]['color'] != 1:
        C += 1
    if G.nodes["Event11"]['color'] != 0:
        C += 1
    if G.nodes["Event12"]['color'] != 3:
        C += 1
    if G.nodes["Event13"]['color'] != 2:
        C += 1
    if G.nodes["Event14"]['color'] != 1:
        C += 1
    if G.nodes["Event15"]['color'] != 0:
        C += 1
    if G.nodes["Event16"]['color'] != 2:
        C += 1
    if G.nodes["Event17"]['color'] != 3:
        C += 1
    if G.nodes["Event22"]['color'] != 3:
        C += 1
    if G.nodes["Event23"]['color'] != 0:
        C += 1
    if G.nodes["Event24"]['color'] != 3:
        C += 1
    if G.nodes["Event25"]['color'] != 1:
        C += 1
    
    #PreferTimes_3
    if G.nodes["Event18"]['color'] != 0:
        C += 1
    #PreferTimes_4
    if G.nodes["Event19"]['color'] != 2:
        C += 1
    #PreferTimes_5
    if G.nodes["Event20"]['color'] != 1:
        C += 1
    #PreferTimes_6
    if G.nodes["Event21"]['color'] != 3:
        C += 1

    return C

def initial_cost_function_den_4pts(G):
    C = 0
    #PreferTimes_3
    if G.nodes["Event18"]['color'] != 0:
        C += 1
    #PreferTimes_4
    if G.nodes["Event19"]['color'] != 2:
        C += 1
    #PreferTimes_5
    if G.nodes["Event20"]['color'] != 1:
        C += 1
    #PreferTimes_6
    if G.nodes["Event21"]['color'] != 3:
        C += 1

    return C

def cost_function_den_25pts(x, G, num_colors):
    coloring = []
    for i in range(len(G)):
        for pos, char in enumerate(x[i*num_colors:(i*num_colors+num_colors)]):
            if int(char):
                coloring.append(pos)

    color_graph_coloring(G, coloring)

    C = 0
    
    if G.nodes["Event1"]['color'] != 0:
        C += 1
    if G.nodes["Event2"]['color'] != 2:
        C += 1
    if G.nodes["Event3"]['color'] != 3:
        C += 1
    if G.nodes["Event4"]['color'] != 1:
        C += 1
    if G.nodes["Event5"]['color'] != 2:
        C += 1
    if G.nodes["Event6"]['color'] != 3:
        C += 1
    if G.nodes["Event7"]['color'] != 3:
        C += 1
    if G.nodes["Event8"]['color'] != 2:
        C += 1
    if G.nodes["Event9"]['color'] != 0:
        C += 1
    if G.nodes["Event10"]['color'] != 1:
        C += 1
    if G.nodes["Event11"]['color'] != 0:
        C += 1
    if G.nodes["Event12"]['color'] != 3:
        C += 1
    if G.nodes["Event13"]['color'] != 2:
        C += 1
    if G.nodes["Event14"]['color'] != 1:
        C += 1
    if G.nodes["Event15"]['color'] != 0:
        C += 1
    if G.nodes["Event16"]['color'] != 2:
        C += 1
    if G.nodes["Event17"]['color'] != 3:
        C += 1
    if G.nodes["Event22"]['color'] != 3:
        C += 1
    if G.nodes["Event23"]['color'] != 0:
        C += 1
    if G.nodes["Event24"]['color'] != 3:
        C += 1
    if G.nodes["Event25"]['color'] != 1:
        C += 1
    
    #PreferTimes_3
    if G.nodes["Event18"]['color'] != 0:
        C += 1
    #PreferTimes_4
    if G.nodes["Event19"]['color'] != 2:
        C += 1
    #PreferTimes_5
    if G.nodes["Event20"]['color'] != 1:
        C += 1
    #PreferTimes_6
    if G.nodes["Event21"]['color'] != 3:
        C += 1

    return C

def cost_function_den_4pts(x, G, num_colors):
    coloring = []
    for i in range(len(G)):
        for pos, char in enumerate(x[i*num_colors:(i*num_colors+num_colors)]):
            if int(char):
                coloring.append(pos)

    color_graph_coloring(G, coloring)

    C = 0
    #PreferTimes_3
    if G.nodes["Event18"]['color'] != 0:
        C += 1
    #PreferTimes_4
    if G.nodes["Event19"]['color'] != 2:
        C += 1
    #PreferTimes_5
    if G.nodes["Event20"]['color'] != 1:
        C += 1
    #PreferTimes_6
    if G.nodes["Event21"]['color'] != 3:
        C += 1

    return C

def cost_function_bra(x, G, classes, teachers, num_colors):
    '''
    DistributeSplitEventsConstraint: At least 1 double lesson(s) Weight=1
    DistributeSplitEventsConstraint: At least 2 double lesson(s) Weight=1
    LimitIdleTimesConstraint: No IDLE times for teachers Weight=3
    ClusterBusyTimesConstraint: Not more than 2 days with lessons Weight=9
    ClusterBusyTimesConstraint: Not more than 3 days with lessons Weight=9
    '''
    pass

def divide_events(events):
    pass

def save_csv(data, nome_csv):
    data_points = pd.DataFrame(data, columns=['Expected Value', 'p', 'Beta0|Gamma|Beta'])
    data_points.to_csv(nome_csv, mode='a', header=False)

    return

def color_graph_num(graph, num_color):
    color_index = 0
    node_list = list(graph.nodes)
    not_allowed_color = []

    # Mark all the vertices as not visited
    visited = {x: False for x in graph.nodes}

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

def color_graph_coloring(graph, coloring):
    for index,node in enumerate(graph.nodes):
        graph.nodes[node]['color'] = coloring[index]

    return

def color_graph_greedy_random(G, alpha):
    
    n = G.number_of_nodes()
    
    #lista de cores
    colors = []
    
    #vertices ordenadas por degree decrescente
    nodes_ordered = sorted(G.degree, key=lambda x: x[1], reverse=True)
    
    
    #lista de cores proibidas
    forbidden_colors = {}
    for (node, d) in nodes_ordered:
        forbidden_colors[node]=[]
        
    color_by_node = {}
    while len(nodes_ordered) > 0:
        node = None
        i = 0
        while i < len(nodes_ordered):
            if random.random() <=  alpha or i == len(nodes_ordered)-1:
                (node, d) = nodes_ordered.pop(i)
                break
            i+=1              
        
        p_colors = list(set.difference(set(colors), set(forbidden_colors[node]))) 
        c = 0
        if len(p_colors) > 0:
            c = p_colors[0]
            color_by_node[node] = c
        else:
            c = len(colors)
            colors.append(c)
            color_by_node[node] = c
        #proibe cor para adjacentes
        for adj in G.neighbors(node):
            forbidden_colors[adj].append(c)

    return color_by_node, colors
 
def create_graph(events):
    G = nx.Graph()
    G.add_nodes_from([(event['Id'], {'color' : None}) for event in events])

    comb = combinations(events, 2)
    for i in comb:
        res0 = set(i[0]['Resources'])
        res1 = i[1]['Resources']
        intersection = [value for value in res0 if value in res1]
        if intersection:
            G.add_edge(i[0]['Id'], i[1]['Id'])
    return G

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

def find_sequence(pos_b, start):
        sequence = []
        end = pos_b[start]
        while end != start:
            sequence.append(end)
            end = pos_b[end]

        sequence.append(start)
        return sequence

def qubits2color(qc, num_nodes, num_colors):
    qdit_ord = []
    color_ord = []
    for i in range(num_nodes):
        for j in range(num_colors):
            pos_a = i*num_colors+j
            pos_b = j*num_nodes+i
            qdit_ord.append(pos_a)
            color_ord.append(pos_b)

    not_visited = set(qdit_ord)

    while not_visited:
        index = next(iter(not_visited))
        sequence = find_sequence(color_ord, index)

        for pos in sequence:
            not_visited.remove(pos)

        if len(sequence) != 1:
            start = sequence.pop()
            while sequence:
                qubit_b = sequence.pop()
                swap(qc[qubit_b], qc[start])
                start = qubit_b

def color2qubits(qc, num_nodes, num_colors):
    qdit_ord = []
    color_ord = []
    for i in range(num_nodes):
        for j in range(num_colors):
            pos_a = i*num_colors+j
            pos_b = j*num_nodes+i
            qdit_ord.append(pos_a)
            color_ord.append(pos_b)

    not_visited = set(qdit_ord)

    while not_visited:
        index = next(iter(not_visited))
        sequence = find_sequence(color_ord, index)

        for pos in sequence:
            not_visited.remove(pos)

        if len(sequence) != 1:
            start = sequence.pop(0)
            while sequence:
                qubit_b = sequence.pop(0)
                swap(qc[qubit_b], qc[start])
                start = qubit_b

def phase_separator(qc, gamma, num_nodes, num_colors):
    qubits2color(qc, num_nodes, num_colors)
    for node in range(num_colors*num_nodes):
        X(qc[node])
    for k in range(num_colors):
        qubits = [k*num_nodes+node for node in range(num_nodes)]
        control = qc[qubits[0]]
        for qub in qubits[1:-1]:
            control = control | qc[qub]
        target = qc[qubits[-1]]
        ctrl(control, RZ, 2*gamma, target)
        #cRz(qc, qubits, 2*gamma)
    for node in range(num_colors*num_nodes):
        X(qc[node])
    color2qubits(qc, num_nodes, num_colors)

def G_gate(p, upper, lower):
    theta = np.arccos(np.sqrt(p))
    ctrl(upper, RY(2*theta), lower)

def dichotomy_tree_gen(num_total):
    dichotomy_tree = []
    n = np.floor(num_total/2)
    m = num_total
    root = n/m
    dichotomy_tree.append((n/m))
    new_values = [(n,m)]

    while new_values:
        leaf = new_values.pop(0)
        if leaf != (0,1) and leaf != (1,1) and leaf != (1,2):
            n = leaf[0]
            m = leaf[1]
            upper_child = (np.floor(n/2), np.floor(m/2))
            lower_child = (np.ceil(n/2), np.ceil(m/2))
            if (upper_child == (1,1) and lower_child == (1,2)) or (lower_child == (1,1) and upper_child == (1,2)):
                temp = upper_child
                upper_child = lower_child
                lower_child = temp
            
            if upper_child == (1,1) or upper_child == (0,1):
                dichotomy_tree.append(None)
            else:
                dichotomy_tree.append(upper_child[0]/upper_child[1])
            if lower_child == (1,1) or lower_child == (0,1):
                dichotomy_tree.append(None)
            else:
                dichotomy_tree.append(lower_child[0]/lower_child[1])
            new_values.append(upper_child)
            new_values.append(lower_child)

    root = build(dichotomy_tree) 
    return root
        
def w_state_preparation(qc):
    n = qc.len()
    X(qc[0])
    if (n & (n-1) == 0) and n != 0:
        for i in range(int(np.log2(n))):
            exp = 2**i
            for j in range(exp):
                ctrl(qc[j], RY, np.pi/2, qc[j+exp])
                cnot(qc[j+exp], qc[j])
    else:
        root = dichotomy_tree_gen(n)
        leafs = []
        G_gate(root.value, qc[0], qc[1])
        cnot(qc[1], qc[0])
        
        target_index = 2
        if root.left.value != None:
            leafs.append((root.left, 0, target_index))
            target_index += 1
        if root.right.value != None:
            leafs.append((root.right, 1, target_index))
            target_index += 1

        while leafs:
            leaf = leafs.pop(0)
            param = leaf[0].value
            upper_qubit_index = leaf[1]
            lower_qubit_index = leaf[2]
            G_gate(param, qc[upper_qubit_index], qc[lower_qubit_index])
            cnot(qc[lower_qubit_index], qc[upper_qubit_index])
            # Left = Upper Child
            upper = leaf[0].left
            if upper != None:
                leafs.append((upper, upper_qubit_index, target_index))
                target_index += 1
            # Right = Lower Child
            lower = leaf[0].right
            if lower != None:
                leafs.append((lower, lower_qubit_index, target_index))
                target_index += 1

def qaoa_min_graph_coloring(p, G, num_colors, beta0, gamma, beta):
    num_nodes = G.number_of_nodes()
    qc = quant((num_nodes*num_colors) + num_nodes)

    # Initial state preparation
    #for i in range(num_nodes):
    #    w_state_preparation(qc[i*num_colors:i*num_colors+num_colors])
    
    coloring = [G.nodes[node]['color'] for node in G.nodes]
    for i, color in enumerate(coloring):
        X(qc[(i*num_colors)+color])

    # Alternate application of operators
    # No need for beta 0 if initial state is W
    mixer(qc, G, beta0, num_nodes, num_colors) # Mixer 0
    for step in range(p):
        phase_separator(qc, gamma[step], num_nodes, num_colors)
        mixer(qc, G, beta[step], num_nodes, num_colors)

    # Measurement
    #result = measure(qc).get()
    return dump(qc)

def qaoa_4pts(par, p, G, num_colors):
    # QAOA parameters
    beta0 = par[0]
    new_par = np.delete(par, 0)
    middle = int(len(par)/2)
    gamma = new_par[:middle]
    beta = new_par[middle:]

    print("Using Following parameters: Beta0:", beta0, "Gamma:", gamma, "Beta:", beta)

    num_nodes = G.number_of_nodes()

    # Dictionary for keeping the results of the simulation
    counts = {}
    # run on local simulator
    result = qaoa_min_graph_coloring(p, G, num_colors, beta0, gamma, beta)
    for i in result.get_states():
        binary = np.binary_repr(i, width=(num_nodes*num_colors)+num_nodes)
        counts[binary] = int(2**20*result.probability(i))

    # Evaluate the data from the simulator
    avr_C       = 0
    min_C       = [0, np.inf]

    for sample in list(counts.keys()):
        if counts[sample] > 0:
            # use sampled bit string x to compute f(x)
            x       = [int(num) for num in list(sample)]
            #tmp_eng = cost_function_den_4pts(x, G, num_colors)
            tmp_eng = cost_function_min(x, G, num_colors)

            # compute the expectation value and energy distribution
            avr_C     = avr_C    + counts[sample]*tmp_eng

            # save best bit string
            if( min_C[1] > tmp_eng):
                min_C[0] = sample
                min_C[1] = tmp_eng

    expectation_value = avr_C/sum(counts.values())

    return expectation_value

def minimization_process_4pts(p, G, num_colors, school):
    data = []
    
    # Initializing QAOA Parameters 
    gamma = [random.uniform(0, 2*np.pi) for _ in range(p)]
    beta0 = random.uniform(0, np.pi)
    beta  = [random.uniform(0, np.pi) for _ in range(p)]

    # Packing and Sending to minimize
    qaoa_par = [beta0]+gamma+beta
    qaoa_args = p, G, num_colors
    print("\nMinimizing function\n")
    res = minimize(qaoa_4pts, qaoa_par, args=qaoa_args, method='Nelder-Mead',
            options={'maxiter': 300, 'disp': True, 'adaptive':True})
    print(res)

    data.append([res['fun'], p, res['x']])
    save_csv(data, "results/"+school+"p"+str(p)+".csv" )

def qaoa_25pts(par, p, G, num_colors):
    # QAOA parameters
    beta0 = par[0]
    new_par = np.delete(par, 0)
    middle = int(len(par)/2)
    gamma = new_par[:middle]
    beta = new_par[middle:]

    print("Using Following parameters: Beta0:", beta0, "Gamma:", gamma, "Beta:", beta)
    
    num_nodes = G.number_of_nodes()

    # Dictionary for keeping the results of the simulation
    counts = {}
    # run on local simulator
    result = qaoa_min_graph_coloring(p, G, num_colors, beta0, gamma, beta)
    for i in result.get_states():
        binary = np.binary_repr(i, width=(num_nodes*num_colors)+num_nodes)
        counts[binary] = int(2**20*result.probability(i))

    # Evaluate the data from the simulator
    avr_C       = 0
    min_C       = [0, np.inf]

    for sample in list(counts.keys()):
        if counts[sample] > 0:
            # use sampled bit string x to compute f(x)
            x       = [int(num) for num in list(sample)]
            tmp_eng = cost_function_den_25pts(x, G, num_colors)
            #tmp_eng = cost_function_min(x, G, num_colors)

            # compute the expectation value and energy distribution
            avr_C     = avr_C    + counts[sample]*tmp_eng

            # save best bit string
            if( min_C[1] > tmp_eng):
                min_C[0] = sample
                min_C[1] = tmp_eng

    expectation_value = avr_C/sum(counts.values())

    return expectation_value

def minimization_process_25pts(p, G, num_colors, school):
    data = []
    
    # Initializing QAOA Parameters 
    gamma = [random.uniform(0, 2*np.pi) for _ in range(p)]
    beta0 = random.uniform(0, np.pi)
    beta  = [random.uniform(0, np.pi) for _ in range(p)]

    # Packing and Sending to minimize
    qaoa_par = [beta0]+gamma+beta
    qaoa_args = p, G, num_colors
    print("\nMinimizing function\n")
    res = minimize(qaoa_25pts, qaoa_par, args=qaoa_args, method='Nelder-Mead',
            options={'maxiter': 300, 'disp': True, 'adaptive':True})
    print(res)

    data.append([res['fun'], p, res['x']])
    save_csv(data, "results/"+school+"p"+str(p)+".csv" )

def create_graphv2(nodes, edges):
    G = nx.Graph()
    G.add_nodes_from([(tuple, {'color' : None}) for tuple in nodes])

    for e, row in enumerate(edges):
        for f, column in enumerate(row):
            if column == 1:
               G.add_edge(nodes[e],nodes[f])

    return G

def cost_function_timetable(x, G, num_colors, list_students):
    coloring = []

    for i in range(len(G)):
        for pos, char in enumerate(x[i*num_colors:(i*num_colors+num_colors)]):
            if int(char):
                coloring.append(pos)

    color_graph_coloring(G, coloring)

    C = 0
    lectures = G.nodes
    for student in list_students:
        lecture_student = [(j,k,l) for (j,k,l) in lectures if l == student]
        timeslots = [G.nodes[key]['color'] for key in list(lecture_student)]
        timeslots.sort()
        new_timeslots = [x%num_colors for x in timeslots]

        # Students shouldn't have lectures at the last time of the day
        for time in new_timeslots:
            if time % num_colors == num_colors-1:
                C += 1

        day = []
        last_lecture = -np.inf
        for time in new_timeslots:
            if last_lecture >= time:
                # Students shouldn't have only one lecture in a day
                if len(day) == 1:
                    C += 1
                #Students shouldn't have many consecutive lectures
                else:
                    for index, lect in enumerate(day):
                        if index+1 < len(day) and day[index+1] == lect+1:
                            C += 1
                day = []

            last_lecture = time
            day.append(last_lecture)
        for index, lect in enumerate(day):
            if index+1 < len(day) and day[index+1] == lect+1:
                C += 1

    return C

def qaoa_first(par, p, G, num_colors, students_list):
    # QAOA parameters
    beta0 = par[0]
    new_par = np.delete(par, 0)
    middle = int(len(par)/2)
    gamma = new_par[:middle]
    beta = new_par[middle:]

    #print("Using Following parameters: Beta0:", beta0, "Gamma:", gamma, "Beta:", beta)

    num_nodes = G.number_of_nodes()

    # Dictionary for keeping the results of the simulation
    counts = {}
    # run on local simulator
    result = qaoa_min_graph_coloring(p, G, num_colors, beta0, gamma, beta)
    for i in result.get_states():
        binary = np.binary_repr(i, width=(num_nodes*num_colors)+num_nodes)
        counts[binary] = int(2**20*result.probability(i))

    # Evaluate the data from the simulator
    avr_C       = 0
    min_C       = [0, np.inf]

    for sample in list(counts.keys()):
        if counts[sample] > 0:
            # use sampled bit string x to compute f(x)
            x       = [int(num) for num in list(sample)]
            tmp_eng = cost_function_timetable(x, G, num_colors, students_list)

            # compute the expectation value and energy distribution
            avr_C     = avr_C    + counts[sample]*tmp_eng

            # save best bit string
            if( min_C[1] > tmp_eng):
                min_C[0] = sample
                min_C[1] = tmp_eng

    expectation_value = avr_C/sum(counts.values())

    return expectation_value

def minimization_process_first(p, G, num_colors, school, students_list):
    # Initializing QAOA Parameters 
    gamma = [random.uniform(0, 2*np.pi) for _ in range(p)]
    beta0 = random.uniform(0, np.pi)
    beta  = [random.uniform(0, np.pi) for _ in range(p)]

    lower_bounds = [0] * ((2*p)+1)
    upper_bounds_beta0 = [np.pi]
    upper_bounds_gamma = [2*np.pi]*p
    upper_bounds_beta  = [np.pi]*p
    upper_bounds = upper_bounds_beta0+upper_bounds_gamma+upper_bounds_beta 

    # Packing and Sending to minimize
    qaoa_par = [beta0]+gamma+beta
    qaoa_args = p, G, num_colors, students_list
    print("\nMinimizing function\n")
    opts = {'bounds' : [lower_bounds, upper_bounds]}#, 'maxiter': 2, } #'maxfevals': 300}
    sigma0 = 2
    
    es = cma.CMAEvolutionStrategy(qaoa_par, sigma0, opts)
    while not es.stop():
        solutions = es.ask()
        es.tell(solutions, [qaoa_first(s, p, G, num_colors, students_list) for s in solutions])
        res = es.result
        print("Saving Results")
        save_csv([[res[1], p, res[0]]], "results/"+school+"p"+str(p)+".csv" )
        #es.disp()
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

def first_example():
    # Problem variables
    num_weeks = 1
    num_days = 1 #5
    num_periods = 6
    num_timeslots = num_days*num_periods

    # Each subject has only one teacher
    # Each teacher teaches only one subject
    num_subjects = 3
    num_teachers = num_subjects
    num_students = 2
    num_rooms = 3

    # Number of features in a room
    # Ex.: has computers, has >40 chairs...
    num_features = 2

    teachers_list = [teacher for teacher in range(num_teachers)]
    students_list = [student for student in range(num_students)]

    #roomFeatures = np.matrix([[0 for feature in range(num_rooms)] for event in range(num_features)])
    roomFeatures = np.matrix([[0, 1, 1],
                             [1, 0, 0]])
    #subjectFeatures = np.matrix([[0 for feature in range(num_features)] for event in range(num_subjects)])
    subjectFeatures = np.matrix([[1, 0],
                                [1, 0],
                                [0, 1]])
    suitableRoom = subjectFeatures*roomFeatures

    # Allocate rooms
    # Each subject will be allocated to the least busy room
    allocations = []
    num_allocations = [0 for room in range(num_rooms)]

    for subject_index, subject in enumerate(suitableRoom.tolist()):
        possible_allocations = []
        for index, room in enumerate(subject):
            if room == 1:
                possible_allocations.append(index)
        print("Subject", subject_index)
        print("Possible Allocations", possible_allocations)

        min_allocations = np.inf
        allocated_room = np.inf
        for alloc_index in possible_allocations:
            if num_allocations[alloc_index] < min_allocations:
                allocated_room = alloc_index
                min_allocations = num_allocations[allocated_room]
        allocations.append((subject_index, allocated_room))
        num_allocations[allocated_room] += 1

    print("\nNumber of Allocations for each Room", num_allocations)
    print("Allocations", allocations)

    # Pair subjects with students
    # lecture = (subject, room, student)
    lectures = [(j,k,l) for j,k in allocations for l in students_list]

    # Generate the lectureConflict Matrix
    # The hard constraints of the problem are included in the matrix
    lectureConflict = [[0 for feature in range(len(lectures))] for event in range(len(lectures))]

    # If two lectures are allocated to the same room,
    # share a student or have the same teacher they
    # cannot be assigned to the same timeslot
    for e, j in enumerate(lectures):
        subject,room,student = j
        for f,a in enumerate(lectures[e+1:]):
            subject2,room2,student2 = a
            if subject == subject2:
                lectureConflict[e][e+1+f] = 1
                lectureConflict[e+1+f][e] = 1
            if student == student2:
                lectureConflict[e][e+1+f] = 1
                lectureConflict[e+1+f][e] = 1
            if room == room2:
                lectureConflict[e][e+1+f] = 1
                lectureConflict[e+1+f][e] = 1

    G = create_graphv2(lectures, lectureConflict)

    return G, students_list

def minimal_example():
    nodes = [('Event1', {'color': None}),
             ('Event2', {'color': None}),
             ('Event3', {'color': None}),
             ('Event4', {'color': None}),
             ('Event5', {'color': None}),
             ('Event6', {'color': None}),
             ('Event7', {'color': None}),
    ]
    edges = [('Event1', 'Event2'),
             ('Event1', 'Event3'),
             ('Event1', 'Event4'),
             ('Event1', 'Event5'),
             ('Event2', 'Event3'),
             ('Event2', 'Event4'),
             ('Event2', 'Event5'),
             ('Event3', 'Event4'),
             ('Event3', 'Event5'),
             ('Event4', 'Event5'),
             ('Event1', 'Event6'),
             ('Event2', 'Event6'),
             ('Event3', 'Event6'),
             ('Event4', 'Event6'),
             ('Event5', 'Event6'),
             ('Event1', 'Event7'),
             ('Event2', 'Event7'),
             ('Event3', 'Event7'),
             ('Event4', 'Event7'),
             ('Event5', 'Event7')
    ]

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    return G

def initial_cost_function_min(G):
    C = 0
    
    if G.nodes["Event1"]['color'] != 1:
        C += 1
    if G.nodes["Event2"]['color'] != 2:
        C += 1
    if G.nodes["Event3"]['color'] != 3:
        C += 1
    if G.nodes["Event4"]['color'] != 4:
        C += 1
    if G.nodes["Event5"]['color'] != 0:
        C += 1
    
    return C

def cost_function_min(x, G, num_colors):
    coloring = []
    for i in range(len(G)):
        for pos, char in enumerate(x[i*num_colors:(i*num_colors+num_colors)]):
            if int(char):
                coloring.append(pos)
    color_graph_coloring(G, coloring)
    
    C = 0
    
    if G.nodes["Event1"]['color'] != 1:
        C += 1
    if G.nodes["Event2"]['color'] != 2:
        C += 1
    if G.nodes["Event3"]['color'] != 3:
        C += 1
    if G.nodes["Event4"]['color'] != 4:
        C += 1
    if G.nodes["Event5"]['color'] != 0:
        C += 1
    
    return C

def main():
    print("Starting program\n")

    # --------------------------
    # School Instances
    # --------------------------
 
    #school = "Den"
    #school = "Bra"
    #school = "Min"
    school = "CEC"

    # Parse XML file
    #events = parseXML('dataset/den-smallschool.xml')
    #events = parseXML('dataset/bra-instance01.xml')
 
    #G = create_graph(events)
    #G = minimal_example()
    G, students_list = first_example()

    # --------------------------
    #  Preparing Conflict Graph
    # --------------------------
    print("\nGraph information")

    #print("Nodes = ", G.nodes)
    #coloring = [G.nodes[node]['color'] for node in G.nodes]
    #print("\nPre-coloring", coloring)

    #degree = [deg for (node, deg) in G.degree()]
    #print("\nDegree of each node", degree)

    # Greedy coloring to be used in cases where a trivial coloring cannot be
    # found
    # ---------------------------------------------------------------------
    #pair = None, G.number_of_nodes(), 0
    #it = 0
    #for i in range (1, 10000):
    #    color_by_node, colors = color_graph_greedy_random(G, 0.7)
    #    if pair[1] > len(colors):
    #        pair = color_by_node, len(colors), it
    #    it+= 1
    # Coloring Graph
    #for key, value in pair[0].items(): 
    #    G.nodes[key]['color'] = value

    #num_colors = pair[1] # Greed coloring colors
    num_colors = 6        # CEC example colors
    #num_colors = 7       # Denmark colors
    #num_colors = 25      # Brazil colors
    print("\nNumber of colors", num_colors)

    # If a suitable coloring can be found without the greedy method use
    # the color_graph_num method
    # -----------------------------------------------------------------
    #color_graph_num(G, num_colors)

    # Minimal example Coloring
    #color_graph_coloring(G, [0,1,2,3,4])
    # First example Coloring
    color_graph_coloring(G, [0,3,1,4,2,5])

    # Verifying Graph consistency
    #----------------------------
    #for i in G.nodes:
    #    print("\nNode",i)
    #    print("Neighbours", G[i])
    #for i in G.nodes:
    #    print("\nNode",i,"Color", G.nodes[i]['color'])
    #    neighbours = [G.nodes[neighbour]['color'] for neighbour in G[i]]
    #    print("Neighbours Colors", neighbours)
    #nx.draw(G, with_labels=True, font_weight='bold')
    #plt.show()

    # Initial Values
    # -------------- 
    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nInitial coloring", coloring)
    #initial_function_value = initial_cost_function_den_25pts(G)
    #print("\nInitial Function Value Max 25:", initial_function_value)
    #initial_function_value = initial_cost_function_den_4pts(G)
    #print("\nInitial Function Value Max 4:", initial_function_value)
    #initial_function_value = initial_cost_function_min(G)
    #print("\nInitial Function Value Max 5:", initial_function_value)

    # -------------
    # Starting QAOA
    # ------------- 
    print("Running QAOA")
    num_nodes = G.number_of_nodes()
    number_of_qubits = num_nodes*num_colors+num_nodes
    print("Necessary number of qubits: ", number_of_qubits)

    # QAOA parameter
    p = 1

    # Minimizing Example CEC
    print("Running minimization process")
    minimization_process_first(p, G, num_colors, school, students_list)
    
    # Minimizing Denmark 4pts
    #print("Running minimization process for 4 points")
    #minimization_process_4pts(p, G, num_colors, school)

    # Minimizing Denmark 25pts
    #print("Running minimization process for 25 points")
    #minimization_process_25pts(p, G, num_colors, school)

if __name__ == '__main__':
    main()