# We import the tools to handle general Graphs
import networkx as nx

import numpy as np
import pprint as pp

# Import miscellaneous tools
from xml_parser import parseXML
from itertools import combinations

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

# Apply the partial mixer for each pair of colors of each node
def mixer(G, beta, num_nodes, num_colors):
    list_nodes = list(G.nodes())
    for u, node in enumerate(G.nodes()):
        print("Verifying Node", node, "in position", u)
        for i in range(num_colors):
            print("Swapping color i", i)
            neighbours_i = neighbourhood(G, num_colors, node, i, list_nodes)
            for j in range(num_colors):
                print("Swapping color j", j)
                if i < j:
                    neighbours_j = neighbourhood(G, num_colors, node, j, list_nodes)
                    neighbours = neighbours_i+neighbours_j

                    #if neighbours == []:
                    #    q_neighbours = None
                    #else:
                    #    q_neighbours = qc[neighbours[0]]
                    #    for neigh in neighbours[1:]:
                    #        q_neighbours = q_neighbours | qc[neigh]
                    #partial_mixer(
                    #        qc,
                    #        q_neighbours,
                    #        qc[num_nodes*num_colors+u],
                    #        qc[i+(num_colors*u)]|qc[j+(num_colors*u)],
                    #        beta)

# Function to check if it is safe to assign color `c` to vertex `v`
def isSafe(graph, color, v, c):
    # check the color of every adjacent vertex of `v`
    for u in graph[v]:
        if color[u] == c:
            return False
 
    return True
 
def kColorable(g, color, k, v, N, states):
 
    # if all colors are assigned, print the solution
    if v == N:
        str1 = " "
        states.append(str1.join([str(color[v]) for v in range(N)]))
        return
 
    # try all possible combinations of available colors
    for c in range(1, k + 1):
        # if it is safe to assign color `c` to vertex `v`
        if isSafe(g, color, v, c):
            # assign color `c` to vertex `v`
            color[v] = c
 
            # recur for the next vertex
            kColorable(g, color, k, v + 1, N, states)
 
            # backtrack
            color[v] = 0

def main():
    print("Starting program\n")
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

    G_tuple = create_graphv2(lectures, lectureConflict)
    G = create_graphv3(lectures, lectureConflict)

    # Graph Information
    print("\nGraph information")

    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nPre-coloring", coloring)

    degree = [deg for (node, deg) in G.degree()]
    print("\nDegree of each node", degree)

    num_colors = num_timeslots
    #num_colors = 6
    print("\nNumber of colors", num_colors)
    initial_coloring = [0,3,1,4,2,5]

    node_list = list(G.nodes)
    #color_graph_num(G, num_colors, node_list[0])
    color_graph_coloring(G, initial_coloring)

    for i in G.nodes:
        print("\nNode",i,"Color", G.nodes[i]['color'])
        neighbours = [G.nodes[neighbour]['color'] for neighbour in G[i]]
        print("Neighbours Colors", neighbours)

    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nInitial coloring", coloring)

    tmp_eng = cost_function_timetable(initial_coloring, G_tuple, num_colors, students_list)
    print("fun", tmp_eng)

    # Probabilities list to be used at the end of the program 
    states = []
    color = [None] * G.number_of_nodes()

    # print all k–colorable configurations of the graph
    kColorable(G, color, num_colors, 0, G.number_of_nodes(), states)

    pp.pprint(states)
    
    print("Running Random Distribution")
    
    # -----------------------------------------------------
    # Evaluate the data from limited number of Measurements
    # -----------------------------------------------------
    measurement_number = 10000
    measurement = np.random.choice(states, measurement_number)
    unique, repet = np.unique(measurement, return_counts=True)
    analysis = dict(zip(unique, repet))

    avr_C       = 0
    hist        = {}
    for sample in list(analysis.keys()):
        # use sampled bit string x to compute f(x)
        x       = [int(num) for num in list(sample.split(" "))]
        tmp_eng = cost_function_timetable(x,G_tuple, num_colors, students_list)

        # compute the expectation value and energy distribution
        avr_C     = avr_C    + analysis[sample]*tmp_eng
        hist[str(round(tmp_eng))] = hist.get(str(round(tmp_eng)),0) + analysis[sample]

    
    print("\nTotal Number of Measurements", measurement_number)
    expected_value = avr_C/measurement_number
    print("Expected Value = ", expected_value)
    print("Objective Function Distribution")
    #pp.pprint(hist)
    #plt.hist(hist)
    #plot_histogram(hist,figsize = (8,6),bar_labels = False)
    #plt.savefig('timetabling-random.pdf') 
    #plt.bar(list(hist.keys()), hist.values(), width=1, color='b')
    #plt.show()

if __name__ == '__main__':
    main()
