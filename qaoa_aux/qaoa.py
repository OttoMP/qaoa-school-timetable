# Import tools for running QAOA
from ket import around, X, H, ctrl, RZ, RX, quant, ket_config, dump, context
from json import load, dump as json_dump

#import math tools
import numpy as np

# We import the tools to handle general Graphs
import networkx as nx
from os import environ

def partial_mixer(neighbour, ancilla, target, beta):
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
                            q_neighbours,
                            qc[num_nodes*num_colors+u],
                            qc[i+(num_colors*u)]|qc[j+(num_colors*u)],
                            beta)

def phase_separator(qc, gamma, num_nodes, num_colors):
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

def qaoa_min_graph_coloring(p, G, num_nodes, num_colors, beta0, gamma, beta, epsilon):
    #ket_config(backend='cpu')
    ket_config(epsilon=epsilon) 
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

def remove_aux_fix_coloring(G, coloring, num_colors):
    # Remove Auxiliar Nodes 
    aux_colors = [4,5]
    for i, color in enumerate(coloring):
        if color == aux_colors[0]:
            coloring[i] = coloring[-2]
        if color == aux_colors[1]:
            coloring[i] = coloring[-1]

    coloring[-2] = aux_colors[0]
    coloring[-1] = aux_colors[1]
    color_graph_from_coloring(G, coloring)

    # Fix Coloring
    for i, node in enumerate(G.nodes):
        if coloring[i] == aux_colors[0] or coloring[i] == aux_colors[1]:
            not_allowed_color = {G.nodes[neighbour]['color'] for neighbour in G[node]}
            if len(not_allowed_color) == num_colors:
                break
            color_index = 0
            while color_index in not_allowed_color:
                color_index = (color_index+1)%num_colors
            coloring[i] = color_index

    coloring[-2] = aux_colors[0]
    coloring[-1] = aux_colors[1]
    color_graph_from_coloring(G, coloring)

def coloring_is_invalid(G):
    #----------------------------
    # Verifying Graph consistency
    #---------------------------
    for e,i in enumerate(G.nodes):
        color = G.nodes[i]['color']
        neighbour_colors = [G.nodes[neighbour]['color'] for neighbour in G[i]]
        if color in neighbour_colors:
            return True
    return False

def qaoa(par, p, initial_G, num_colors, epsilon, cost_function, school, it_num):
    ctx = context()
    # --------------------------
    # Unpacking QAOA parameters
    # --------------------------
    beta0 = par[0]
    middle = int(len(par)/2)
    gamma = par[1:middle+1]
    beta = par[middle+1:]
    num_nodes = initial_G.number_of_nodes()

    # --------------------------
    # Verifying Parameters
    # --------------------------
    #print("Using Following parameters: Beta0:", beta0, "Gamma:", gamma, "Beta:", beta, "Epsilon:", epsilon)

    # --------------------------
    # Running QAOA on simulator
    # --------------------------
    G = nx.Graph()
    G.add_nodes_from(initial_G)
    G.add_edges_from(initial_G.edges)
    initial_coloring = [initial_G.nodes[node]['color'] for node in initial_G.nodes]
    color_graph_from_coloring(G, initial_coloring)
    
    result = qaoa_min_graph_coloring(p, initial_G, num_nodes, num_colors, beta0, gamma, beta, epsilon)
    
    #print("Number of States", len(result.get_states()))
    #print("State Vector", result.show('b6:b6:b6:b6:b6:b6'))

    # --------------------------
    # Counting resulting states
    # --------------------------
    counts = {} # Dictionary for keeping the results of the simulation
    for i in result.states:
        binary = f'{i:0{(num_nodes*num_colors)+num_nodes}b}'
        counts[binary] = int(2**20*result.probability(i))
    
    # --------------------------
    # Evaluate the data from the simulator
    # --------------------------
    avr_function_value = 0
    valid_solutions = {}
    invalid_solutions = {}
    for sample in list(counts.keys()):
        if counts[sample] > 0:
            # extracting x to compute f(x)
            # ----------------------------
            x       = [int(num) for num in list(sample)]
            
            # Coloring Graph with counts[sample]
            # ----------------------------------
            coloring = []
            for i in range(len(G)):
                for pos, char in enumerate(x[i*num_colors:(i*num_colors+num_colors)]):
                    if int(char):
                        coloring.append(pos)
            color_graph_from_coloring(G, coloring)

            # Removing auxiliary nodes for correct function evaluation
            # --------------------------------------------------------
            remove_aux_fix_coloring(G, coloring, num_colors)
            if coloring_is_invalid(G):
                invalid_solutions[sample] = counts[sample]
                continue
            else:
                valid_solutions[sample] = counts[sample]

            # Computing fx
            # ------------
            fx = cost_function(G)

            # Compute the expectation value and energy distribution
            # -----------------------------------------------------
            avr_function_value = avr_function_value + counts[sample]*fx

    # Return expected value
    expected_valid = avr_function_value/sum(valid_solutions.values())
    with open(f'timing/{school}_{p}_{it_num}.txt', 'a') as file:
        exec_time = ctx.get_return('exec_time')
        file.write(f'{exec_time}, ')

    return expected_valid

def color_graph_from_coloring(graph, coloring):
    for index,node in enumerate(graph.nodes):
        graph.nodes[node]['color'] = coloring[index]
    return
