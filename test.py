# We import the tools to handle general Graphs
import networkx as nx

import numpy as np
import pprint as pp

# Import miscellaneous tools
from xml_parser import parseXML
from itertools import combinations

# Import tools for running QAOA
import random
import cma
from scipy.optimize import minimize
from ket import *
from ket.lib import swap, within
from ket.gates.more import u3

# We import plotting tools
import pandas as pd
import matplotlib.pyplot as plt
from qiskit.visualization import plot_histogram

def create_graphv1(events):
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

def create_graphv2(nodes, edges):
    G = nx.Graph()
    G.add_nodes_from([(tuple, {'color' : None}) for tuple in nodes])

    for e, row in enumerate(edges):
        for f, column in enumerate(row):
            if column == 1:
               G.add_edge(nodes[e],nodes[f])

    return G

def create_graphv3(nodes, edges):
    G = nx.Graph()
    G.add_nodes_from([(num, {'color' : None}) for num in range(len(nodes))])

    for e, row in enumerate(edges):
        for f, column in enumerate(row):
            if column == 1:
                G.add_edge(e,f)

    return G

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

def initial_cost_function_timetable(x, G, num_colors, list_students):

    color_graph_coloring(G, x)

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
    print("Neighbourhood function")
    neighbours = list(G[node])
    neighbours_index = [list_nodes.index(neigh) for neigh in neighbours]
    print("Neighbours from node", node,":", neighbours)

    neighbours_color_qubit = [color+(num_colors*u) for u in neighbours_index]
    print("Neighbour qubits", neighbours_color_qubit)

    return neighbours_color_qubit

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

def qaoa_first(par, p, G, num_colors, students_list):
    # QAOA parameters
    beta0 = par[0]
    new_par = np.delete(par, 0)
    middle = int(len(new_par)/2)
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
    data = []
    
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
    print(lower_bounds)
    print(upper_bounds)
    opts = {'maxiter': 2,'bounds' : [lower_bounds, upper_bounds], } # 'maxfevals': 300}
    sigma0 = 1
    #es = cma.CMAEvolutionStrategy(qaoa_par, sigma0, opts)
    #es.optimize(qaoa_first, args=qaoa_args)
    #res = es.result

    res = cma.fmin(qaoa_first, qaoa_par, sigma0, args=qaoa_args, options=opts)
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

def main():

    school = "CEC"
    
    G, students_list = first_example()
    
    num_colors = 6 #First example colors
    print("\nNumber of colors", num_colors)
    
    color_graph_coloring(G, [0,3,1,4,2,5])
    
    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nInitial coloring", coloring)

    init_fun = initial_cost_function_timetable(coloring, G, num_colors, students_list)
    print("Initial Function Value", init_fun)

    # -------------
    # Starting QAOA
    # ------------- 
    print("Running QAOA")
    number_of_qubits = G.number_of_nodes()*num_colors+G.number_of_nodes()
    print("Necessary number of qubits: ", number_of_qubits)

    # QAOA parameter
    p = 1

    # Minimizing Example CEC
    print("Running minimization process")
    minimization_process_first(p, G, num_colors, school, students_list)
    
if __name__ == '__main__':
    main()