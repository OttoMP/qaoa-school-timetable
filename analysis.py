# Import tools for running QAOA
from ket import *

#import math tools
import numpy as np

# We import the tools to handle general Graphs
import networkx as nx

# Import miscellaneous tools
import random, sys
from xml_parser import parseXML
from itertools import combinations
import pandas as pd
import pprint as pp

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

def cost_function_den_25pts(G):
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

def cost_function_den_4pts(G):
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

def color_graph_greedy(G):
    pair = None, G.number_of_nodes(), 0
    it = 0
    for i in range (1, 10000):
        color_by_node, colors = color_graph_greedy_aux(G, 0.7)
        if pair[1] > len(colors):
            pair = color_by_node, len(colors), it
        it+= 1
    # Coloring Graph
    for key, value in pair[0].items(): 
        G.nodes[key]['color'] = value

def color_graph_greedy_aux(G, alpha):
    
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
    G.add_nodes_from([(event['Id'], {'color' : None}) for event in events])

    comb = combinations(events, 2)
    for i in comb:
        res0 = set(i[0]['Resources'])
        res1 = i[1]['Resources']
        intersection = [value for value in res0 if value in res1]
        if intersection:
            G.add_edge(i[0]['Id'], i[1]['Id'])
    return G

def add_auxiliary_nodes(G, num_colors):
    #Adding Auxiliary nodes for mixing 
    nodes = [('Aux1', {'color': num_colors}),
             ('Aux2', {'color': num_colors+1}),
    ]
    G.add_nodes_from(nodes)
    
    for n in G.nodes:
        if n != 'Aux1' and n != 'Aux2':
            G.add_edge(n, 'Aux1')
            G.add_edge(n, 'Aux2')

def instance_den4():
    # --------------------------
    # Parse XML file
    # --------------------------
    events = parseXML('dataset/den-smallschool.xml')

    # --------------------------
    #  Preparing Conflict Graph
    # --------------------------
    initial_G = create_graph_from_events(events)
    # --------------------------
    #  Coloring Conflict Graph
    # --------------------------
    #color_graph_greedy(initial_G)

    # Coloring 23 points 
    initial_coloring =  [1, 0, 2, 3, 1, 2, 1, 2, 3, 0, 0, 2, 0, 3, 1, 3, 0, 1, 0, 3, 2, 2, 1, 2, 3]
    color_graph_from_coloring(initial_G, initial_coloring)
    
    initial_function_value = cost_function_den_25pts(initial_G)
    print("\nInitial Function Value Max 25", initial_function_value)
    initial_function_value = cost_function_den_4pts(initial_G)
    print("\nInitial Function Value Max 4", initial_function_value)

    # Adding extras nodes to allow mixing
    add_auxiliary_nodes(initial_G, 4)

    return initial_G

def instance_den5():
    # --------------------------
    # Parse XML file
    # --------------------------
    events = parseXML('dataset/den-smallschool.xml')

    # --------------------------
    #  Preparing Conflict Graph
    # --------------------------
    initial_G = create_graph_from_events(events)
    # --------------------------
    #  Coloring Conflict Graph
    # --------------------------
    
    # Greedy coloring to be used in cases where a trivial coloring cannot be
    # found
    # -----------------------------------------------------------------
    #color_graph_greedy(initial_G)
    
    # If a suitable coloring can be found without the greedy method use
    # the color_graph_num method
    # -----------------------------------------------------------------
    num_colors = 5 
    color_graph_from_num(initial_G, num_colors)
    #initial_coloring =  [1, 0, 2, 3, 1, 2, 1, 2, 3, 0, 0, 2, 0, 3, 1, 3, 0, 1, 0, 3, 2, 2, 1, 2, 3]
    #color_graph_from_coloring(initial_G, initial_coloring)
    
    initial_function_value = cost_function_den_25pts(initial_G)
    print("\nInitial Function Value Max 25", initial_function_value)
    initial_function_value = cost_function_den_4pts(initial_G)
    print("\nInitial Function Value Max 4", initial_function_value)

    return initial_G

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

def show_results_manual(school, p, initial_G, num_colors):
    #----------------------------
    # Loading QAOA parameters 
    #----------------------------
    print("\nMin Expected Value: ", 3.5)
    
    beta0 = 1.17814704
    gamma = [4.87014693]
    beta =  [2.3562754]
    print("Using Following parameters:")
    print("Beta0:", beta0)
    print("Gamma:", gamma)
    print("Beta:", beta)

    # -------------
    # Starting QAOA
    # ------------- 
    print("\nRunning QAOA")
    num_nodes = initial_G.number_of_nodes()
    number_of_qubits = num_nodes*num_colors+num_nodes
    print("Necessary number of qubits: ", number_of_qubits)
    # --------------------------
    # Running QAOA on simulator
    # --------------------------
    G = nx.Graph()
    G.add_nodes_from(initial_G)
    G.add_edges_from(initial_G.edges)
    
    result = qaoa_min_graph_coloring(p, initial_G, num_nodes, num_colors, beta0, gamma, beta)
    
    # --------------------------
    # Counting resulting states
    # --------------------------
    counts = {} # Dictionary for keeping the results of the simulation
    states = []
    probabilities = []
    for i in result.states:
        binary = f'{i:0{(num_nodes*num_colors)+num_nodes}b}'
        counts[binary] = int((2**20)*result.probability(i))
        states.append(binary)
        probabilities.append(result.probability(i))

    print("Number of States", len(counts))
    #pp.pprint(counts)

    # ------------------------- 
    # Evaluate the data from the simulator
    # --------------------------
    print('\n --- SIMULATION RESULTS ---')
    avr_function_value = 0
    min_function_value = [0, np.inf]
    hist        = {}

    for sample in list(counts.keys()):
        if counts[sample] > 0:
            # use sampled bit string x to compute f(x)
            x = [int(num) for num in list(sample)]
            
            # Coloring Graph with counts[sample]
            coloring = []
            for i in range(len(G)):
                for pos, char in enumerate(x[i*num_colors:(i*num_colors+num_colors)]):
                    if int(char):
                        coloring.append(pos)
            color_graph_from_coloring(G, coloring)
            
            if  school == "Den4":
                remove_aux_fix_coloring(G, coloring, num_colors)
                if coloring_is_invalid(G):
                    print("Invalid coloring", coloring)
                    print("Probability", counts[sample])
                #fx = cost_function_den_25pts(G)
                fx = cost_function_den_4pts(G)
            elif school == "Den5":
                #fx = cost_function_den_25pts(G)
                fx = cost_function_den_4pts(G)

            # compute the expectation value and energy distribution
            avr_function_value = avr_function_value + counts[sample]*fx
            hist[str(round(fx))] = hist.get(str(round(fx)),0) + counts[sample]

            # save best bit string
            if( min_function_value[1] > fx):
                min_function_value[0] = sample
                min_function_value[1] = fx

    print('\n--- Function Distribution  ---')
    total_counts = sum(counts.values())
    print("Total Number of Measurements", total_counts)
    expected_value = avr_function_value/total_counts
    print("Expected Value = ", expected_value)
    print("Objective Function Distribution")
    pp.pprint(hist)

    # -----------------------------------------------------
    # Evaluate the data from limited number of Measurements
    # -----------------------------------------------------
    measurement_number = 10000
    measurement = np.random.choice(states, measurement_number, p=probabilities)
    unique, repet = np.unique(measurement, return_counts=True)
    analysis = dict(zip(unique, repet))

    avr_function_value       = 0
    hist        = {}
    for sample in list(analysis.keys()):
        # use sampled bit string x to compute f(x)
        x = [int(num) for num in list(sample)]
        
        # Coloring Graph with counts[sample]
        coloring = []
        for i in range(len(G)):
            for pos, char in enumerate(x[i*num_colors:(i*num_colors+num_colors)]):
                if int(char):
                    coloring.append(pos)
        color_graph_from_coloring(G, coloring)

        if coloring_is_invalid(G):
            print("Invalid coloring", coloring)
            print("Probability", counts[sample])

        if  school == "Den4":
            remove_aux_fix_coloring(G, coloring, num_colors)
            #fx = cost_function_den_25pts(G)
            fx = cost_function_den_4pts(G)
        elif school == "Den5":
            #fx = cost_function_den_25pts(G)
            fx = cost_function_den_4pts(G)

        # compute the expectation value and energy distribution
        avr_function_value = avr_function_value + analysis[sample]*fx
        hist[str(round(fx))] = hist.get(str(round(fx)),0) + analysis[sample]

    print('\n--- Function Distribution  ---')
    print("Total Number of Measurements", measurement_number)
    expected_value = avr_function_value/measurement_number
    print("Expected Value = ", expected_value)
    print("Objective Function Distribution")
    pp.pprint(hist)
    
    # -----------------------------------------------------
    print('--------------------------')
    print('--- Individual States  ---\n')
    print("Best result found: ", min_function_value[0])
    print("Number of times result showed: ", counts[min_function_value[0]])
    print("Percentage of times result showed: ", (counts[min_function_value[0]]/total_counts)*100)
    print("Objective function value: ", min_function_value[1])

    list_qubits = min_function_value[0]
    best_coloring = []
    for i in range(len(G)):
        for pos, char in enumerate(list_qubits[i*num_colors:(i*num_colors+num_colors)]):
            if int(char):
                # color = pos
                best_coloring.append(pos)

    remove_aux_fix_coloring(G, best_coloring, num_colors)
    print("Best Coloring",best_coloring)
    #print("Best Coloring Qudits values")
    #for i in range(len(G)):
    #    print(list_qubits[i*num_colors:(i*num_colors+num_colors)])
    print('\n')
    
    #-----------------------------
    max_counts = max(counts, key=lambda key: counts[key])
    print("Most commom result found: ", max_counts)
    print("Number of times result showed: ", counts[max_counts])
    print("Percentage of times result showed: ", (counts[max_counts]/total_counts)*100)
    # Coloring Graph with max_counts
    coloring = []
    for i in range(len(G)):
        for pos, char in enumerate(max_counts[i*num_colors:(i*num_colors+num_colors)]):
            if int(char):
                coloring.append(pos)
    
    color_graph_from_coloring(G, coloring)

    if  school == "Den4":
        remove_aux_fix_coloring(G, best_coloring, num_colors)
        #common_value = cost_function_den_25pts(G)
        common_value = cost_function_den_4pts(G)
    elif school == "Den5":
        #common_value = cost_function_den_25pts(G)
        common_value = cost_function_den_4pts(G)
    
    print("Objective function value: ", common_value)

    list_qubits = max_counts
    common_coloring = []
    for i in range(len(G)):
        for pos, char in enumerate(list_qubits[i*num_colors:(i*num_colors+num_colors)]):
            if int(char):
                # color = pos
                common_coloring.append(pos)

    print("Most Common Coloring",common_coloring)
    #print("\nMost Common Coloring Qudits values")
    #for i in range(len(G)):
    #    print(list_qubits[i*num_colors:(i*num_colors+num_colors)])

def print_graph_info(initial_G):
    print("--------------------------")
    print("Graph information\n")
    print("Nodes = ", initial_G.nodes)
    degree = [deg for (node, deg) in initial_G.degree()]
    print("\nDegree of each node", degree)
    initial_coloring = [initial_G.nodes[node]['color'] for node in initial_G.nodes]
    print("\nInitial coloring", initial_coloring)
    num_colors = len(set(initial_coloring))
    print("\nNumber of colors", num_colors)
    
    #----------------------------
    # Verifying Graph consistency
    #---------------------------
    print("----------------------------")
    print("Verifying Graph consistency")
    for e,i in enumerate(initial_G.nodes):
        print("\nNode",e,"Color", initial_G.nodes[i]['color'])
        color_and_neighbour = [(neighbour, initial_G.nodes[neighbour]['color']) for neighbour in initial_G[i]]
        print("Neighbours | Color")
        for pair in color_and_neighbour:
            print(pair)

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
            if len(not_allowed_color) >= num_colors:
                continue
            color_index = 0
            while color_index in not_allowed_color:
                color_index = (color_index+1)%num_colors
            coloring[i] = color_index

    coloring[-2] = aux_colors[0]
    coloring[-1] = aux_colors[1]
    color_graph_from_coloring(G, coloring)


def main():
    print("Starting program\n")

    # QAOA parameter
    p = int(sys.argv[1])

    school = "Den4"
    #school = "Den5"
    print("Analysing instance: ", school)

    if  school == "Den4":
        initial_G = instance_den4()
        coloring = [initial_G.nodes[node]['color'] for node in initial_G.nodes]
        num_colors = len(set(coloring))
    elif school == "Den5":
        initial_G = instance_den5()
        coloring = [initial_G.nodes[node]['color'] for node in initial_G.nodes]
        #num_colors = len(set(coloring))
        num_colors = 5
    
    print_graph_info(initial_G) 

    num_nodes = initial_G.number_of_nodes()
    number_of_qubits = num_nodes*num_colors+num_nodes
    print("Number of Nodes", num_nodes)
    print("Number of colors", num_colors)
    print("Necessary number of qubits: ", number_of_qubits)

    show_results_manual(school, p, initial_G, num_colors)
    #show_results_cma()
    #show_results_cobyla() 

if __name__ == '__main__':
    main()