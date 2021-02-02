#import math tools
import numpy as np

# We import the tools to handle general Graphs
import networkx as nx

# Import miscellaneous tools
import pprint as pp
import time
import progressbar

# We import plotting tools
import pandas as pd
import matplotlib.pyplot as plt
from   matplotlib import cm
from   matplotlib.ticker import LinearLocator, FormatStrFormatter

# Import tools for running QAOA
import random
from scipy.optimize import minimize, fmin, Bounds
from sympy.combinatorics.graycode import GrayCode
from ket.lib import swap

def show_figure(fig):
    new_fig = plt.figure()
    new_mngr = new_fig.canvas.manager
    new_mngr.canvas.figure = fig
    fig.set_canvas(new_mngr.canvas)
    plt.show(fig)

# Compute the value of the cost function
def cost_function_C(x,G,num_colors):
    C = 0
    for color in range(num_colors):
        color_used = 0
        for node in G.nodes():
            if x[node*num_colors + color] == 1:
                color_used = 1
                break
        C += color_used

    return C

def save_csv(data, nome_csv):
    data_points = pd.DataFrame(data, columns=['Expected Value', 'p', 'Graph Number'])
    data_points.to_csv(nome_csv, mode='a', header=False)

    return

def create_graph_line():
    # Line
    # {"0": [1], "1": [0,2], "2": [1,3], "3": [2]}
    n     = 4
    V     = np.arange(0,n,1)
    E     =[(0,1),(1,2),(2,3)] 

    G     = nx.Graph()
    G.add_nodes_from([(vertex, {'color' : None}) for vertex in V])
    G.add_edges_from(E)

    # Generate plot of the Graph
    colors       = ['r' for node in G.nodes()]
    default_axes = plt.axes(frameon=True)
    pos          = nx.spring_layout(G)

    return G

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
                                if (graph.nodes[neighbour]['color'] != None 
                                    and graph.nodes[neighbour]['color'] < num_color) ]
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

def create_graph(nodes, edges):
    G = nx.Graph()
    G.add_nodes_from([(num, {'color' : None}) for num in range(len(nodes))])

    for e, row in enumerate(edges):
        for f, column in enumerate(row):
            if column == 1:
                G.add_edge(e,f)

    return G

def create_graph_tuple(nodes, edges):
    G = nx.Graph()
    G.add_nodes_from([(tuple, {'color' : None}) for tuple in nodes])

    for e, row in enumerate(edges):
        for f, column in enumerate(row):
            if column == 1:
               G.add_edge(nodes[e],nodes[f])

    return G

def cRz(qc, qubits, gamma):
        num_qubits = len(qubits)
        qubits_t = qubits[1:]
        exp = np.power(2, num_qubits-1)
        a = GrayCode(num_qubits)
        gray_list = list(a.generate_gray())[1:]

        # Add the necessary gates following the Gray Code
        u1(gamma/exp, qc[qubits[0]])
        ctrl(qc[qubits[0]], x, qc[qubits[1]])
        counter = 1
        for i, node in enumerate(qubits_t):
            for j in range(np.power(2, i+1)-1):
                u1(np.power(-1, counter)*gamma/exp, qc[node])
                counter += 1
                codes = zip(gray_list[counter-1][::-1], gray_list[counter][::-1])
                enumerator = 0
                for u,v in codes:
                    if u != v:
                        ctrl(qc[qubits[enumerator]], x, qc[node])
                    enumerator += 1
            if i < len(qubits_t)-1:
                u1(np.power(-1, counter)*gamma/exp, qc[node])
                counter += 1
                ctrl(qc[node], x, qc[qubits_t[i+1]])
        u1(np.power(-1, counter)*gamma/exp, qc[qubits_t[-1]])

def toffoli(qc, controls, target):
        all_qubits = controls+[target]
        h(qc[target])
        cRz(qc, all_qubits, np.pi)
        h(qc[target])

def partial_mixer(qc, neighbour, ancilla, target, beta):
        if neighbour == []:
            x(qc[ancilla])
        else:
            for node in neighbour:
                x(qc[node])
            toffoli(qc, neighbour, ancilla)
            for node in neighbour:
                x(qc[node])

        ## Phase correction
        #u1(-beta, qc[ancilla])

        # Controlled Rxx
        h(qc[target[0]])
        h(qc[target[1]])
        ctrl(qc[target[0]], x, qc[target[1]])
        cRz(qc, [ancilla]+[target[1]], 2*beta)
        ctrl(qc[target[0]], x, qc[target[1]])
        h(qc[target[0]])
        h(qc[target[1]])

        # Controlled Ryy
        rx(-np.pi/2,qc[target[0]])
        rx(-np.pi/2,qc[target[1]])
        ctrl(qc[target[0]], x, qc[target[1]])
        cRz(qc, [ancilla]+[target[1]], 2*beta)
        ctrl(qc[target[0]], x, qc[target[1]])
        rx(np.pi/2,qc[target[0]])
        rx(np.pi/2,qc[target[1]])

        if neighbour == []:
            x(qc[ancilla])
        else:
            for node in neighbour:
                x(qc[node])
            toffoli(qc, neighbour, ancilla)
            for node in neighbour:
                x(qc[node])

def neighbourhood(G, num_colors, node, color):
    neighbour = G[node]
    neighbour_qubit = [color+(num_colors*u) for u in neighbour]

    return neighbour_qubit

# Apply the partial mixer for each pair of colors of each node
def mixer(qc, G, beta, num_nodes, num_colors):
    for u in range(num_nodes):
        for i in range(num_colors):
            neighbours_i = neighbourhood(G, num_colors, u, i)
            for j in range(num_colors):
                if i < j:
                    neighbours_j = neighbourhood(G, num_colors, u, j)
                    neighbours = neighbours_i+neighbours_j
                    partial_mixer(
                            qc,
                            neighbours,
                            num_nodes*num_colors+u,
                            [i+(num_colors*u), j+(num_colors*u)],
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
        x(qc[node])
    for k in range(num_colors):
        qubits = [k*num_nodes+node for node in range(num_nodes)]
        cRz(qc, qubits, 2*gamma)
    for node in range(num_colors*num_nodes):
        x(qc[node])
    color2qubits(qc, num_nodes, num_colors)

def qaoa_min_graph_coloring(p, G, num_colors, gamma, beta0, beta):
    num_nodes = G.number_of_nodes()
    qc = quant((num_nodes*num_colors) + num_nodes)

    # Initial state preparation
    coloring = [G.nodes[node]['color'] for node in G.nodes]
    for i, color in enumerate(coloring):
        x(qc[(i*num_colors)+color])

    # Alternate application of operators
    #mixer(qc, G, beta0, num_nodes, num_colors) # Mixer 0
    for step in range(p):
        #phase_separator(qc, G, gamma[step], num_nodes, num_colors)
        mixer(qc, G, beta[step], num_nodes, num_colors)
    #print(report())
    # Measurement
    result = measure(qc).get()
    #result = dump(qc)
    return result

def main():
    print("Starting program")
    G = create_graph_line()

    # Graph Information
    print("\nGraph information")

    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nPre-coloring", coloring)

    degree = [deg for (node, deg) in G.degree()]
    print("\nDegree of each node", degree)

    #num_colors = max(degree)+2
    num_colors = 4
    print("\nNumber of colors", num_colors)

    node_list = list(G.nodes)
    color_graph_num(G, num_colors, node_list[0])

    for i in G.nodes:
        print("\nNode",i,"Color", G.nodes[i]['color'])
        neighbours = [G.nodes[neighbour]['color'] for neighbour in G[i]]
        print("Neighbours Colors", neighbours)

    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nInitial coloring", coloring)

    print("Running QAOA")
    p = 1
    gamma = [random.uniform(0, 2*np.pi) for _ in range(p)]
    beta0 =  random.uniform(0, np.pi)
    beta  = [random.uniform(0, np.pi) for _ in range(p)]
    num_nodes = G.number_of_nodes()

    print("Using Following parameters:")
    print("Beta 0:", beta0)
    print("Gamma:", gamma)
    print("Beta:", beta)
    print("\n")

    result = qaoa_min_graph_coloring(p, G, num_colors, gamma, beta0, beta)
    print("==Result==")
    #print(result.show())
    shots = 100
    counts = {}

    # run on local simulator
    for _ in progressbar.progressbar(range(shots)):
        result = qaoa_min_graph_coloring(p, G, num_colors, gamma, beta0, beta)
        binary = np.binary_repr(result, width=(num_nodes*num_colors)+num_nodes)
        if binary in counts:
            counts[binary] += 1
        else:
            counts[binary] = 1

    pp.pprint(counts)

    # Evaluate the data from the simulator
    avr_C       = 0
    min_C       = [0, G.number_of_nodes()+1]

    for sample in list(counts.keys()):
        # use sampled bit string x to compute C(x)
        x         = [int(num) for num in list(sample)]
        tmp_eng   = cost_function_C(x,G,num_colors)

        # compute the expectation value and energy distribution
        avr_C     = avr_C    + counts[sample]*tmp_eng

        # save best bit string
        if( min_C[1] > tmp_eng):
            min_C[0] = sample
            min_C[1] = tmp_eng

    M1_sampled   = avr_C/shots
    print('\n --- SIMULATION RESULTS ---\n')
    final_coloring = []
    list_qubits = min_C[0]
    print(list_qubits)

    for i in range(len(G)):
        for pos, char in enumerate(list_qubits[i*num_colors:(i*num_colors+num_colors)]):
            if int(char):
                # color = pos
                final_coloring.append(pos)

    print("\nFinal Coloring",final_coloring)
    print("\nFinal Coloring Qudits values")
    for i in range(len(G)):
        print(list_qubits[i*num_colors:(i*num_colors+num_colors)])

    #print('\nThe approximate solution is x* = %s with C(x*) = %d' % (min_C[0],min_C[1]))
    #print('The number of times this solution showed was: %d \n' %(counts[min_C[0]]))
    #print('The sampled mean value is Mp_sampled = %.02f' % (M1_sampled))

    print("New Graph information")
    print("\nDegree of each node", degree)
    print("\nNumber of colors", num_colors)
    color_graph_coloring(G, final_coloring)
    for i in G.nodes:
        print("\nNode",i,"Color", G.nodes[i]['color'])
        neighbours = [G.nodes[neighbour]['color'] for neighbour in G[i]]
        print("Neighbours Colors", neighbours)

    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nFinal coloring", final_coloring)

if __name__ == '__main__':
    main()