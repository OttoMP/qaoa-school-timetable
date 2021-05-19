#import math tools
import numpy as np

# We import the tools to handle general Graphs
import networkx as nx

# Import miscellaneous tools
from xml_parser import parseXML
from itertools import combinations
import pprint as pp
import progressbar
from binarytree import build

# We import plotting tools
import pandas as pd
import matplotlib.pyplot as plt
from qiskit.visualization import plot_histogram

# Import tools for running QAOA
import random
from scipy.optimize import minimize
from ket import *
from ket.lib import swap, within
from ket.gates.more import u3

# Parallelization tools
import ray
import multiprocessing

#num_cores = multiprocessing.cpu_count()
num_cores = 2
ray.init(num_cpus=num_cores)

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
        print("Testing Node", u, "as root")
        success = True
        queue.append(u)
        visited[u] = True

        while queue:
            # Dequeue a vertex from queue and color it
            source = queue.pop(0)
            not_allowed_color = {graph.nodes[neighbour]['color'] for neighbour in graph[source]
                                    if (graph.nodes[neighbour]['color'] != None) }
            print("Source", source)
            print("Not allowed color", not_allowed_color)
            if len(not_allowed_color) == num_color:
                print("Node", u, "not correct\n")
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
            print("Succesfully Colored\n")
            break

    return

def color_graph_coloring(graph, coloring):
    for index,node in enumerate(graph.nodes):
        graph.nodes[node]['color'] = coloring[index]

    return

def color_graph_greedy(G):
    node_list = list(G.nodes)
    V = G.number_of_nodes()
    result = [-1] * V
 
    # Assign the first color to first vertex
    root = node_list[0]
    G.nodes[root]['color'] = 0
 
    # A temporary array to store the available colors.
    # True value of available[cr] would mean that the
    # color cr is assigned to one of its adjacent vertices
    available = [False for _ in range(V)]
 
    # Assign colors to remaining V-1 vertices
    for u in node_list[1:]:
         
        # Process all adjacent vertices and
        # flag their colors as unavailable
        not_allowed_color = [G.nodes[neighbour]['color'] for neighbour in G[u]
                                if (G.nodes[neighbour]['color'] != None) ]
        for color in not_allowed_color:
            available[color] = True
 
        # Find the first available color
        cr = 0
        while cr < V:
            if (available[cr] == False):
                break
            cr += 1
             
        # Assign the found color
        G.nodes[u]['color'] = cr
 
        # Reset the values back to false
        # for the next iteration
        available = [False for _ in range(V)]

def create_graphv2(nodes, edges):
    G = nx.Graph()
    G.add_nodes_from([(num, {'color' : None}) for num in nodes])

    for e, row in enumerate(edges):
        for f, column in enumerate(row):
            if column == 1:
                G.add_edge(nodes[e],nodes[f])

    return G

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
    RY(theta, lower)
    cnot(upper, lower)
    RY(-theta, lower)
    cnot(upper, lower)

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

def qaoa(par, p, G, num_colors):
    # QAOA parameters
    beta0 = par.pop(0)
    middle = int(len(par)/2)
    gamma = par[:middle]
    beta = par[middle:]

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

@ray.remote
def minimization_process(p, G, num_colors, school):
    data = []
    
    # Initializing QAOA Parameters 
    gamma = [random.uniform(0, 2*np.pi) for _ in range(p)]
    beta0 = random.uniform(0, np.pi) for _ in range(p)
    beta  = [random.uniform(0, np.pi) for _ in range(p)]

    # Packing and Sending to minimize
    qaoa_par = [beta0]+gamma+beta
    qaoa_args = p, G, num_colors
    print("\nMinimizing function\n")
    res = minimize(qaoa, qaoa_par, args=qaoa_args, method='Nelder-Mead',
            options={'maxiter': 1, 'maxfev': 1, 'disp': True, 'adaptive':True})
    #print(res)

    data.append([res['fun'], p, res['x']])
    save_csv(data, "results/"+school+"p"+str(p)+".csv" )

    return [['fun'], p, res['x']]

def main():
    print("Starting program\n")

    # ----------------------------------
    # Minimal Example Preparation Begins
    # ----------------------------------
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
        #print("Subject", subject_index)
        #print("Possible Allocations", possible_allocations)

        min_allocations = np.inf
        allocated_room = np.inf
        for alloc_index in possible_allocations:
            if num_allocations[alloc_index] < min_allocations:
                allocated_room = alloc_index
                min_allocations = num_allocations[allocated_room]
        allocations.append((subject_index, allocated_room))
        num_allocations[allocated_room] += 1

    #print("\nNumber of Allocations for each Room", num_allocations)
    #print("Allocations", allocations)

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

    # ----------------------------------
    # Minimal Example Preparation Ends
    # ----------------------------------

    # School Instances
    # Parse XML file
    events = parseXML('dataset/den-smallschool.xml')
    #events = parseXML('dataset/bra-instance01.xml')

    school = "Den"
    #school = "Bra"

    #G = create_graphv2(lectures, lectureConflict)
    G = create_graph(events)

    # --------------------------
    #  Preparing Conflict Graph
    # --------------------------
    print("\nGraph information")

    print("Nodes = ", G.nodes)
    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nPre-coloring", coloring)

    degree = [deg for (node, deg) in G.degree()]
    print("\nDegree of each node", degree)

    num_colors = 5
    print("\nNumber of colors", num_colors)

    color_graph_num(G, num_colors)
    #color_graph_coloring(G, initial_coloring)
    #color_graph_greedy(G)

    for i in G.nodes:
        print("\nNode",i,"Color", G.nodes[i]['color'])
        neighbours = [G.nodes[neighbour]['color'] for neighbour in G[i]]
        print("Neighbours Colors", neighbours)

    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nInitial coloring", coloring)

    #nx.draw(G, with_labels=True, font_weight='bold')
    #plt.show()

    print("Running QAOA")
    number_of_qubits = G.number_of_nodes()*num_colors+G.number_of_nodes()
    print("Necessary number of qubits: ", number_of_qubits)

    # QAOA parameter
    p = 1

    # Parallel task using ray
    expected_value_sample = ray.get([minimization_process.remote(p, G, num_colors, school) for iteration in progressbar.progressbar(range(1))])

if __name__ == '__main__':
    main()