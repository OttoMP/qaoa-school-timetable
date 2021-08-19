# Import tools for running QAOA
import random
from ket import *

#import math tools
import numpy as np

# We import the tools to handle general Graphs
import networkx as nx

# Import miscellaneous tools
import gc
import os
import psutil
from pympler import asizeof, muppy, summary, tracker, refbrowser
import sys
from xml_parser import parseXML
from itertools import combinations
import pprint as pp

def cost_function_timetable(G, num_colors, list_students):
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

def phase_separator_ad_hoc(qc, gamma, num_nodes, num_colors):
    #if G.nodes["Event1"]['color'] != 0:
    RZ(2*gamma, qc[(num_colors*0)+0])
    #if G.nodes["Event2"]['color'] != 2:
    RZ(2*gamma, qc[(num_colors*1)+2])
    #if G.nodes["Event3"]['color'] != 3:
    RZ(2*gamma, qc[(num_colors*2)+3])
    #if G.nodes["Event4"]['color'] != 1:
    RZ(2*gamma, qc[(num_colors*3)+1])
    #if G.nodes["Event5"]['color'] != 2:
    RZ(2*gamma, qc[(num_colors*4)+2])
    #if G.nodes["Event6"]['color'] != 3:
    RZ(2*gamma, qc[(num_colors*5)+3])
    #if G.nodes["Event7"]['color'] != 3:
    RZ(2*gamma, qc[(num_colors*6)+3])
    #if G.nodes["Event8"]['color'] != 2:
    RZ(2*gamma, qc[(num_colors*7)+2])
    #if G.nodes["Event9"]['color'] != 0:
    RZ(2*gamma, qc[(num_colors*8)+0])
    #if G.nodes["Event10"]['color'] != 1:
    RZ(2*gamma, qc[(num_colors*9)+1])
    #if G.nodes["Event11"]['color'] != 0:
    RZ(2*gamma, qc[(num_colors*10)+0])
    #if G.nodes["Event12"]['color'] != 3:
    RZ(2*gamma, qc[(num_colors*11)+3])
    #if G.nodes["Event13"]['color'] != 2:
    RZ(2*gamma, qc[(num_colors*12)+2])
    #if G.nodes["Event14"]['color'] != 1:
    RZ(2*gamma, qc[(num_colors*13)+1])
    #if G.nodes["Event15"]['color'] != 0:
    RZ(2*gamma, qc[(num_colors*14)+0])
    #if G.nodes["Event16"]['color'] != 2:
    RZ(2*gamma, qc[(num_colors*15)+2])
    #if G.nodes["Event17"]['color'] != 3:
    RZ(2*gamma, qc[(num_colors*16)+3])
    #if G.nodes["Event22"]['color'] != 3:
    RZ(2*gamma, qc[(num_colors*21)+3])
    #if G.nodes["Event23"]['color'] != 0:
    RZ(2*gamma, qc[(num_colors*22)+0])
    #if G.nodes["Event24"]['color'] != 3:
    RZ(2*gamma, qc[(num_colors*23)+3])
    #if G.nodes["Event25"]['color'] != 1:
    RZ(2*gamma, qc[(num_colors*24)+1])
    
    #PreferTimes_3
    #if G.nodes["Event18"]['color'] != 0:
    RZ(2*gamma, qc[(num_colors*17)+0])
    #PreferTimes_4
    #if G.nodes["Event19"]['color'] != 2:
    RZ(2*gamma, qc[(num_colors*18)+2])
    #PreferTimes_5
    #if G.nodes["Event20"]['color'] != 1:
    RZ(2*gamma, qc[(num_colors*19)+1])
    #PreferTimes_6
    #if G.nodes["Event21"]['color'] != 3:
    RZ(2*gamma, qc[(num_colors*20)+3])

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

def qaoa_min_graph_coloring(p, G, num_nodes, num_colors, beta0, gamma, beta, tr):
    #print("Memory Usage start result", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff() 
    # --------------------------
    # Initializing qubits
    # --------------------------
    qc = quant((num_nodes*num_colors) + num_nodes)
    #print("Memory Usage after initi qubit", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff() 

    # --------------------------
    # Initial state preparation
    # --------------------------
    coloring = [G.nodes[node]['color'] for node in G.nodes]
    for i, color in enumerate(coloring):
        X(qc[(i*num_colors)+color])
    #print("Memory Usage after state prep", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff() 

    # --------------------------
    # Alternate application of operators
    # --------------------------
    mixer(qc, G, beta0, num_nodes, num_colors) # Mixer 0
    for step in range(p):
        phase_separator(qc, gamma[step], num_nodes, num_colors)
        #phase_separator_ad_hoc(qc, gamma[step], num_nodes, num_colors)
        mixer(qc, G, beta[step], num_nodes, num_colors)
    #print("Memory Usage after operators", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff() 

    # --------------------------
    # Measurement
    # --------------------------
    #result = measure(qc).get()
    return dump(qc)

def qaoa(par, p, G, num_colors, students_list, tr):
    #print("Memory Usage start qaoa call", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff() 
    # --------------------------
    # Unpacking QAOA parameters
    # --------------------------
    beta0 = par[0]
    middle = int(len(par)/2)
    gamma = par[1:middle+1]
    beta = par[middle+1:]

    #print("Memory Usage after parameters", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff() 
    # --------------------------
    # Verifying Parameters
    # --------------------------
    #print("Using Following parameters:\nBeta0:", beta0, "\nGamma:", gamma, "\nBeta:", beta)

    num_nodes = G.number_of_nodes()
    # --------------------------
    # Running QAOA on simulator
    # --------------------------
    #print("Memory Usage before result", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff() 
    result = qaoa_min_graph_coloring(p, G, num_nodes, num_colors, beta0, gamma, beta, tr)
    #print("Memory Usage after result", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff() 
    
    #print("Number of States", len(result.get_states()))
    #print("State Vector", result.show('b6:b6:b6:b6:b6:b6'))

    # --------------------------
    # Counting resulting states
    # --------------------------
    counts = {} # Dictionary for keeping the results of the simulation
    for i in result.states:
        binary = f'{i:0{(num_nodes*num_colors)+num_nodes}b}'
        counts[binary] = int(2**20*result.probability(i))
    #counts = counting_states(result, num_nodes, num_colors)
    #print("Memory Usage after get states", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff() 
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

            fx = cost_function_timetable(G, num_colors, students_list)
            #fx = cost_function_den_25pts(G)
            #fx = cost_function_den_4pts(G)

            # compute the expectation value and energy distribution
            avr_function_value = avr_function_value + counts[sample]*fx

            # save best bit string
            if( min_function_value[1] > fx):
                min_function_value[0] = sample
                min_function_value[1] = fx
    
    #print("Memory Usage after counts", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff() 

    expected_value = avr_function_value/sum(counts.values())
    # Return expected value
    #print("Memory Usage end qaoa call", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff()
    
    return expected_value

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

def create_graph_from_list(nodes, edges):
    G = nx.Graph()
    G.add_nodes_from([(tuple, {'color' : None}) for tuple in nodes])

    for e, row in enumerate(edges):
        for f, column in enumerate(row):
            if column == 1:
               G.add_edge(nodes[e],nodes[f])

    return G

def create_full_graph(tr):
    #print("Memory Usage start graph", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff() 
    # --------------------------
    # Parse XML file
    # --------------------------
    events = parseXML('dataset/den-smallschool.xml')
    #print("Memory Usage after parse", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff() 

    # --------------------------
    #  Preparing Conflict Graph
    # --------------------------
    G = create_graph_from_events(events)
    #print("Memory Usage after creating graph", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff() 
    
    #print("--------------------------")
    #print("Graph information\n")
    
    #print("Nodes = ", G.nodes)
    #coloring = [G.nodes[node]['color'] for node in G.nodes]
    #print("\nPre-coloring", coloring)

    #degree = [deg for (node, deg) in G.degree()]
    #print("\nDegree of each node", degree)

    # --------------------------
    #  Coloring Conflict Graph
    # --------------------------
    
    # Greedy coloring to be used in cases where a trivial coloring cannot be
    # found
    # -----------------------------------------------------------------
    #color_graph_greedy(G)
    
    # If a suitable coloring can be found without the greedy method use
    # the color_graph_num method
    # -----------------------------------------------------------------
    num_colors = 5 # Denmark colors
    color_graph_from_num(G, num_colors)
    #print("Memory Usage after coloring graph", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    #tr.print_diff() 
    
    #coloring = [G.nodes[node]['color'] for node in G.nodes]
    #coloring =  [1, 0, 2, 3, 1, 2, 1, 2, 3, 0, 0, 2, 0, 3, 1, 3, 0, 1, 0, 3, 2, 2, 1, 2, 3]
    #coloring =  [0, 2, 3, 1, 2, 3, 3, 2, 0, 1, 0, 3, 2, 1, 0, 2, 3, 0, 2, 1, 3, 3, 0, 3, 1]
    #color_graph_from_coloring(G, coloring)
    
    #num_colors = len(set(coloring))
    #print("\nNumber of colors", num_colors)
    #print("\nInitial coloring", coloring)

    #initial_function_value = cost_function_den_25pts(G)
    #print("\nInitial Function Value Max 25", initial_function_value)
    #initial_function_value = cost_function_den_4pts(G)
    #print("\nInitial Function Value Max 4", initial_function_value)

    # ---------------------------
    # Verifying Graph consistency
    #----------------------------
    #print("----------------------------")
    #print("Verifying Graph consistency")
    #for i in G.nodes:
        #print("\nNode",i,"Color", G.nodes[i]['color'])
        #color_and_neighbour = [(neighbour, G.nodes[neighbour]['color']) for neighbour in G[i]]
        #print("Neighbours | Color")
        #for pair in color_and_neighbour:
            #print(pair)
    return G, num_colors
 
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

        min_allocations = np.inf
        allocated_room = np.inf
        for alloc_index in possible_allocations:
            if num_allocations[alloc_index] < min_allocations:
                allocated_room = alloc_index
                min_allocations = num_allocations[allocated_room]
        allocations.append((subject_index, allocated_room))
        num_allocations[allocated_room] += 1

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

    G = create_graph_from_list(lectures, lectureConflict)

    return G, students_list

def main():
    print("Starting program\n")
    tr = tracker.SummaryTracker()
    tr.print_diff() 
    # --------------------------
    # School Instances
    # --------------------------
    school = "CEC"
           
    #----------------------------
    # Starting QAOA
    #----------------------------
    initial_G, students_list = first_example()
    #initial_G, num_colors = create_full_graph(tr)
    
    num_colors = 6        # CEC example colors
    color_graph_from_coloring(initial_G, [0,3,1,4,2,5])
    G = nx.Graph()
    G.add_nodes_from(initial_G)
    G.add_edges_from(initial_G.edges)
    initial_coloring = [initial_G.nodes[node]['color'] for node in initial_G.nodes]
    color_graph_from_coloring(G, initial_coloring)
    
    num_nodes = G.number_of_nodes()
    number_of_qubits = num_nodes*num_colors+num_nodes
    print("Necessary number of qubits: ", number_of_qubits)

    print("----------------------------")
    print("Running QAOA")
    # QAOA parameter
    p = int(sys.argv[1])

    gamma = [random.uniform(0, 2*np.pi) for _ in range(p)]
    beta0 = random.uniform(0, np.pi)
    beta  = [random.uniform(0, np.pi) for _ in range(p)]
    s = [beta0]+gamma+beta
    tr.print_diff() 
    for i in range(10):
        print("Iteration:", i)
        print("Memory Usage before qaoa call", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
        print("Muppy all objects", len(muppy.get_objects()))
        tr.print_diff()
        print("Garbage counter", len(gc.garbage))
        print()
        #ev = qaoa(s, p, G, num_colors, tr)
        ev = qaoa(s, p, G, num_colors, students_list, tr)
        print("Memory Usage after qaoa call", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
        print("Muppy all objects", len(muppy.get_objects()))
        tr.print_diff() 
        print("Garbage counter", len(gc.garbage))
        print("-------------------")
    print("Memory Usage after gc", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    print("Muppy all objects", len(muppy.get_objects()))
    tr.print_diff() 
    print("Garbage counter", len(gc.get_objects()))

if __name__ == '__main__':
    main()