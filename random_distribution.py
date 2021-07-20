# We import the tools to handle general Graphs
import networkx as nx

# Import miscellaneous tools
import numpy as np
import pprint as pp
from xml_parser import parseXML
from itertools import combinations

# We import plotting tools
import matplotlib.pyplot as plt

def cost_function_timetable(x, G, num_colors, list_students):

    color_graph_from_coloring(G, x)

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

def color_graph_from_coloring(graph, coloring):
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

def create_graph_from_list(nodes, edges):
    G = nx.Graph()
    G.add_nodes_from([(tuple, {'color' : None}) for tuple in nodes])

    for e, row in enumerate(edges):
        for f, column in enumerate(row):
            if column == 1:
               G.add_edge(nodes[e],nodes[f])

    return G

def first_example():
    print("--------------------------")
    print("Starting graph preparation\n")
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
        print("Possible Room Allocations", possible_allocations)

        min_allocations = np.inf
        allocated_room = np.inf
        for alloc_index in possible_allocations:
            if num_allocations[alloc_index] < min_allocations:
                allocated_room = alloc_index
                min_allocations = num_allocations[allocated_room]
        allocations.append((subject_index, allocated_room))
        num_allocations[allocated_room] += 1

    print("\nNumber of Allocations for each Room", num_allocations)
    print("Final Allocations (subject, room)", allocations)

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

    print("\nCreating Graph\n")
    G_tuple = create_graph_from_list(lectures, lectureConflict)
    G = create_graph(lectures, lectureConflict)

    return G, G_tuple, students_list

def main():
    print("Starting program\n")

    # --------------------------
    # School Instances
    # --------------------------
    school = "CEC"

    # --------------------------
    # Parse XML file
    # --------------------------
    #events = parseXML('dataset/den-smallschool.xml')

    # --------------------------
    #  Preparing Conflict Graph
    # --------------------------
    G, G_tuple, students_list = first_example()
    
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
    num_colors = 6        # CEC example colors
    print("\nNumber of colors", num_colors)
    
    color_graph_from_coloring(G, [0,3,1,4,2,5])
    initial_coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nInitial coloring", initial_coloring)

    # ---------------------------
    # Verifying Graph consistency
    #----------------------------
    print("----------------------------")
    print("Verifying Graph consistency")
    for e,i in enumerate(G.nodes):
        print("\nNode",e,"Color", G.nodes[i]['color'])
        color_and_neighbour = [(neighbour, G.nodes[neighbour]['color']) for neighbour in G[i]]
        print("Neighbours | Color")
        for pair in color_and_neighbour:
            print(pair)
            
    tmp_eng = cost_function_timetable(initial_coloring, G_tuple, num_colors, students_list)
    print("fun", tmp_eng)

    # Probabilities list to be used at the end of the program 
    states = []
    color = [None] * G.number_of_nodes()

    # print all k–colorable configurations of the graph
    kColorable(G, color, num_colors, 0, G.number_of_nodes(), states)

    #pp.pprint(states)
    
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
    pp.pprint(hist)
    #plt.hist(hist)
    #plot_histogram(hist,figsize = (8,6),bar_labels = False)
    #plt.savefig('timetabling-random.pdf') 
    #plt.bar(list(hist.keys()), hist.values(), width=1, color='b')
    #plt.show()

if __name__ == '__main__':
    main()
