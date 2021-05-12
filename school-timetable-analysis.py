#import math tools
import numpy as np

# We import the tools to handle general Graphs
import networkx as nx

# Import miscellaneous tools
from xml_parser import parseXML
from itertools import combinations
import pprint as pp
import progressbar

# We import plotting tools
import pandas as pd
import matplotlib.pyplot as plt
from qiskit.visualization import plot_histogram

# Import tools for running QAOA
import random
from scipy.optimize import minimize
from ket import *
from ket.plugins import make_quantum
from ket.lib import swap, within

# Parallelization tools
import multiprocessing

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

def cost_function_den(x, G, num_colors):
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
    if G.nodes["Event19"]['color'] != 1:
        C += 1
    #PreferTimes_5
    if G.nodes["Event20"]['color'] != 2:
        C += 1
    #PreferTimes_6
    if G.nodes["Event21"]['color'] != 3:
        C += 1

    return C

def save_csv(data, nome_csv):
    data_points = pd.DataFrame(data, columns=['Expected Value', 'p', 'Beta0|Gamma|Beta'])
    data_points.to_csv(nome_csv, mode='a', header=False)

    return

def color_graph_num(graph, num_color, root):
    color_index = 0
    not_allowed_color = []

    # Mark all the vertices as not visited
    visited = {x: False for x in graph.nodes}

    # Create a queue for BFS
    queue = []

    # Mark the source node as
    # visited and enqueue it
    queue.append(root)
    visited[root] = True

    while queue:
        # Dequeue a vertex from queue and color it
        source = queue.pop(0)
        not_allowed_color = [graph.nodes[neighbour]['color'] for neighbour in graph[source]
                                if (graph.nodes[neighbour]['color'] != None) ]
        while color_index in not_allowed_color:
            color_index = (color_index+1)%num_color
        graph.nodes[source]['color'] = color_index
        not_allowed_color = []

        # Get all adjacent vertices of the
        # dequeued vertex s. If a adjacent
        # has not been visited, then mark it
        # visited and enqueue it
        for i in graph[source]:
            if visited[i] == False:
                queue.append(i)
                visited[i] = True

    return

def color_graph_coloring(graph, coloring):
    for index,node in enumerate(graph.nodes):
        graph.nodes[node]['color'] = coloring[index]

    return

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

def neighbourhood(G, num_colors, node, color):
    neighbour = G[node]
    neighbour_qubit = [color+(num_colors*u) for u, neigh in enumerate(neighbour)]

    return neighbour_qubit

# Apply the partial mixer for each pair of colors of each node
def mixer(qc, G, beta, num_nodes, num_colors):
    for u, node in enumerate(G.nodes):
        for i in range(num_colors):
            neighbours_i = neighbourhood(G, num_colors, node, i)
            for j in range(num_colors):
                if i < j:
                    neighbours_j = neighbourhood(G, num_colors, node, j)
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

def phase_separator(qc, G, gamma, num_nodes, num_colors):
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

@make_quantum
def G_gate(p, target):
    matrix = [[np.sqrt(p), -np.sqrt(1-p)],[np.sqrt(1-p), np.sqrt(p)]]
    return matrix*target

def dichotomy_tree_gen(num_total):
    dichotomy_tree = []
    n = np.floor(num_total/2)
    m = num_total
    root = n/m
    dichotomy_tree.append(n/m)
    new_values = [(n,m)]

    while new_values:
        leaf = new_values.pop(0)
        n = leaf[0]
        m = leaf[1]
        upper_child = (np.floor(n/2), np.floor(m/2))
        lower_child = (np.ceil(n/2), np.ceil(m/2))
        if upper_child != (0,1) and upper_child != (1,1):
            new_values.append(upper_child)
            dichotomy_tree.append(upper_child[0]/upper_child[1])
        if lower_child != (0,1) and lower_child != (1,1):
            new_values.append(lower_child)
            dichotomy_tree.append(lower_child[0]/lower_child[1])
    
    return dichotomy_tree
        
        
def w_state_preparation(qc, num_colors, num_nodes):
    n = qc.len()
    X(qc[0])
    if (n & (n-1) == 0) and n != 0:
        for i in range(int(np.log2(n))):
            exp = 2**i
            for j in range(exp):
                ctrl(qc[j], RY, np.pi/2, qc[j+exp])
                cnot(qc[j+exp], qc[j])
    else:
        dich_tree = dichotomy_tree_gen(n)
        for i in range(int(ceil(np.log2(n)))):
            exp = 2**i
            for j in range(exp):
                if j+exp > n:
                    break
                param = dich_tree.pop(0)
                ctrl(qc[j], G_gate, param, qc[j+exp])
                cnot(qc[j+exp], qc[j])

def w4_state_preparation(qc, num_colors, num_nodes):
    for qudit in range(num_nodes):
        X(qc[qudit*num_colors])
        ctrl(qc[qudit*num_colors], RY, np.pi/2, qc[qudit*num_colors+1])
        cnot(qc[qudit*num_colors+1], qc[qudit*num_colors])
        ctrl(qc[qudit*num_colors], RY, np.pi/2, qc[qudit*num_colors+2])
        cnot(qc[qudit*num_colors+2], qc[qudit*num_colors])
        ctrl(qc[qudit*num_colors+1], RY, np.pi/2, qc[qudit*num_colors+3])
        cnot(qc[qudit*num_colors+3], qc[qudit*num_colors+1])

def qaoa_min_graph_coloring(p, G, num_colors, gamma, beta0, beta):
    num_nodes = G.number_of_nodes()
    qc = quant((num_nodes*num_colors) + num_nodes)

    # Initial state preparation
    w4_state_preparation(qc, num_colors, num_nodes)
    
    #coloring = [G.nodes[node]['color'] for node in G.nodes]
    #for i, color in enumerate(coloring):
    #    X(qc[(i*num_colors)+color])

    # Alternate application of operators
    mixer(qc, G, beta0, num_nodes, num_colors) # Mixer 0
    for step in range(p):
        phase_separator(qc, G, gamma[step], num_nodes, num_colors)
        mixer(qc, G, beta[step], num_nodes, num_colors)

    # Measurement
    #result = measure(qc).get()
    return dump(qc)


def qaoa(par, p, G, num_colors):
    # QAOA parameters
    beta0, par= par[0], par[1:]
    middle = int(len(par)/2)
    gamma = par[:middle]
    beta = par[middle:]

    num_nodes = G.number_of_nodes()

    # Dictionary for keeping the results of the simulation
    counts = {}
    # run on local simulator
    result = qaoa_min_graph_coloring(p, G, num_colors, gamma, beta0, beta)
    for i in result.get_states():
        binary = np.binary_repr(i, width=(num_nodes*num_colors)+num_nodes)
        counts[binary] = int(2**20*result.probability(i))

    # Evaluate the data from the simulator
    avr_C       = 0
    min_C       = [0, 9999]

    for sample in list(counts.keys()):
        if counts[sample] > 0:
            # use sampled bit string x to compute f(x)
            x       = [int(num) for num in list(sample)]
            #tmp_eng = cost_function_timetable(x, G, num_colors, students_list)
            tmp_eng = cost_function_den(x, G, num_colors)

            # compute the expectation value and energy distribution
            avr_C     = avr_C    + counts[sample]*tmp_eng

            # save best bit string
            if( min_C[1] > tmp_eng):
                min_C[0] = sample
                min_C[1] = tmp_eng

    expectation_value = avr_C/sum(counts.values())

    return expectation_value

def minimization_process(p, G, num_colors, school):
    data = []
    gamma = [random.uniform(0, 2*np.pi) for _ in range(p)]
    beta0 = [random.uniform(0, np.pi)]
    beta  = [random.uniform(0, np.pi) for _ in range(p)]
    qaoa_par = beta0+gamma+beta
    qaoa_args = p, G, num_colors
    print("\nMinimizing function\n")
    res = minimize(qaoa, qaoa_par, args=qaoa_args, method='Nelder-Mead',
            options={'maxiter': 300, 'disp': True, 'adaptive':True})
    #print(res)

    data.append([res['fun'], p, res['x']])
    save_csv(data, "results/"+school+"p"+str(p)+".csv" )

    return [['fun'], p, res['x']]

def main():
    print("Starting program\n")

    # QAOA parameter
    p = 1

    school = "Den"
    #school = "Bra"

    # Reading values from CSV file
    result_list = pd.read_csv("results/"+school+"p"+str(p)+".csv", header=0)

    # Parsing necessary values
    expected_value_list = list(result_list["Expected Value"])
    qaoa_par_list = list(result_list["Beta0|Gamma|Beta"])
    min_expected_value = min(expected_value_list)
    min_expected_value_index = expected_value_list.index(min_expected_value)
    min_qaoa_par = qaoa_par_list[min_expected_value_index][1:-1].split()
    
    print("Expected Value List")
    print(expected_value_list)
    print("Min Expected Value")
    print(min_expected_value, '\n')
    print("Beta0|Gamma|Beta")
    print(qaoa_par_list, '\n')

    beta0 = float(min_qaoa_par[0])
    gamma = [float(par) for par in min_qaoa_par[1:p+1]]
    beta  = [float(par) for par in min_qaoa_par[p+1:]]

    print("Using Following parameters:")
    print("Beta 0:", beta0)
    print("Gamma:", gamma)
    print("Beta:", beta)
    print("\n")

    # Preparing Conflict Graph for QAOA Analysis 
    # parse xml file
    events = parseXML('dataset/den-smallschool.xml')
    #events = parseXML('dataset/bra-instance01.xml')
    G = create_graph(events)

    # Graph Information
    #print("\nGraph information")

    coloring = [G.nodes[node]['color'] for node in G.nodes]
    #print("\nPre-coloring", coloring)

    degree = [deg for (node, deg) in G.degree()]
    #print("\nDegree of each node", degree)

    num_colors = 5
    #print("\nNumber of colors", num_colors)

    node_list = list(G.nodes)
    #color_graph_num(G, num_colors, node_list[0])
    #color_graph_coloring(G, initial_coloring)

    for i in G.nodes:
    #    print("\nNode",i,"Color", G.nodes[i]['color'])
        neighbours = [G.nodes[neighbour]['color'] for neighbour in G[i]]
    #    print("Neighbours Colors", neighbours)

    coloring = [G.nodes[node]['color'] for node in G.nodes]
    #print("\nInitial coloring", coloring)

    #nx.draw(G, with_labels=True, font_weight='bold')
    #plt.show()

    # Starting QAOA
    print("Running QAOA")
    number_of_qubits = G.number_of_nodes()*num_colors+G.number_of_nodes()
    print("Necessary number of qubits: ", number_of_qubits)
    
    num_nodes = G.number_of_nodes()
    # Dictionary for keeping the results of the simulation
    counts = {}
    # run on local simulator
    #result = qaoa_min_graph_coloring(p, G, num_colors, gamma, beta0, beta)
    #for i in result.get_states():
    #    binary = np.binary_repr(i, width=(num_nodes*num_colors)+num_nodes)
    #    prob = int(2**20*result.probability(i))
    #    if prob > 0:
    #        counts[binary] = prob
    #pp.pprint(counts)

    #qc = quant(4)
    #qc = quant(1)
    #X(qc[0])
    #G_gate(0.5, qc[0])
    ##w_state_preparation(qc, 4, 1)
    #result = dump(qc)
    #for i in result.get_states():
        #binary = np.binary_repr(i, width=4)
        #prob = int(2**20*result.probability(i))
        #if prob > 0:
            #counts[binary] = prob
    #pp.pprint(counts)

    dich_tree = dichotomy_tree_gen(11)
    print("Tree", dich_tree)
    dich_tree = dichotomy_tree_gen(6)
    print("Tree", dich_tree)
    '''
    print("==Result==")

    # Evaluate the data from the simulator
    avr_C       = 0
    min_C       = [0, 9999]
    hist        = {}

    for k in range(num_colors+1):
        hist[str(k)] = hist.get(str(k),0)

    for sample in list(counts.keys()):
        if counts[sample] > 0:
            # use sampled bit string x to compute f(x)
            x         = [int(num) for num in list(sample)]
            #tmp_eng   = cost_function_timetable(x, G, num_colors, students_list)
            tmp_eng = cost_function_den(x, G, num_colors)

            # compute the expectation value and energy distribution
            avr_C     = avr_C    + counts[sample]*tmp_eng
            hist[str(round(tmp_eng))] = hist.get(str(round(tmp_eng)),0) + counts[sample]

            # save best bit string
            if( min_C[1] > tmp_eng):
                min_C[0] = sample
                min_C[1] = tmp_eng

    total_counts = sum(counts.values())
    print("Total Number of Measurements", total_counts)
    expected_value = avr_C/total_counts
    print("Expected Value = ", expected_value)

    print('\n --- SIMULATION RESULTS ---\n')
    print("Best result found: ", min_C[0])
    print("Number of times result showed: ", counts[min_C[0]])
    print("Percentage of times result showed: ", (counts[min_C[0]]/total_counts)*100)
    print("Objective function value: ", min_C[1])
    print()

    list_qubits = min_C[0]
    best_coloring = []
    for i in range(len(G)):
        for pos, char in enumerate(list_qubits[i*num_colors:(i*num_colors+num_colors)]):
            if int(char):
                # color = pos
                best_coloring.append(pos)

    print("\nBest Coloring",best_coloring)
    print("\nBest Coloring Qudits values")
    for i in range(len(G)):
        print(list_qubits[i*num_colors:(i*num_colors+num_colors)])

    print("\nNew Graph information")
    print("\nDegree of each node", degree)
    #print("\nNumber of colors", num_colors)
    color_graph_coloring(G, best_coloring)
    for i in G.nodes:
        print("\nNode",i,"Color", G.nodes[i]['color'])
        neighbours = [G.nodes[neighbour]['color'] for neighbour in G[i]]
        print("Neighbours Colors", neighbours)

    print('\n')
    #-----------------------------
    max_counts = max(counts, key=lambda key: counts[key])
    print("Most commom result found: ", max_counts)
    print("Number of times result showed: ", counts[max_counts])
    print("Percentage of times result showed: ", (counts[max_counts]/total_counts)*100)
    #max_value = cost_function_timetable(max_counts, G, num_colors, students_list)
    max_value = cost_function_den(x, G, num_colors)
    print("Objective function value: ", max_value)

    list_qubits = max_counts
    maximum_coloring = []
    for i in range(len(G)):
        for pos, char in enumerate(list_qubits[i*num_colors:(i*num_colors+num_colors)]):
            if int(char):
                # color = pos
                maximum_coloring.append(pos)

    print("\nMost Common Coloring",maximum_coloring)
    print("\nMost Common Coloring Qudits values")
    for i in range(len(G)):
        print(list_qubits[i*num_colors:(i*num_colors+num_colors)])

    print("\nNew Graph information")
    print("\nDegree of each node", degree)
    #print("\nNumber of colors", num_colors)
    color_graph_coloring(G, maximum_coloring)
    for i in G.nodes:
        print("\nNode",i,"Color", G.nodes[i]['color'])
        neighbours = [G.nodes[neighbour]['color'] for neighbour in G[i]]
        print("Neighbours Colors", neighbours)

    #-----------------------------
    '''
    '''
    print("Histogram", hist)
    hist_max = sum(counts.values())
    hist_3 = (hist['3']/hist_max)*100
    hist_2 = (hist['2']/hist_max)*100
    hist_0 = (hist['0']/hist_max)*100
    hist_1 = (hist['1']/hist_max)*100
    hist_4 = (hist['4']/hist_max)*100
    hist_5 = (hist['5']/hist_max)*100
    print("Histogram 3", hist_3)
    print("Histogram 2", hist_2)
    print("Histogram 1", hist_1)
    print("Histogram 0", hist_0)
    print("Histogram 4", hist_4)
    print("Histogram 5", hist_5)
    #print('The cost function is distributed as: \n')
    #plot_histogram(hist,figsize = (8,6),bar_labels = False)
    #plt.savefig("histogram.pdf")
    '''
if __name__ == '__main__':
    main()