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
    if G.nodes["Event19"]['color'] != 1:
        C += 1
    #PreferTimes_5
    if G.nodes["Event20"]['color'] != 2:
        C += 1
    #PreferTimes_6
    if G.nodes["Event21"]['color'] != 3:
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

def create_graph_from_events(events):
    G = nx.Graph()
    event_list = [event['Id'] for event in events]
    #G.add_nodes_from([(event['Id'], {'color' : None}) for event in events])
    G.add_nodes_from([(num, {'color' : None}) for num,event in enumerate(events)])

    comb = combinations(events, 2)
    for i in comb:
        res0 = set(i[0]['Resources'])
        res1 = i[1]['Resources']
        intersection = [value for value in res0 if value in res1]
        if intersection:
            G.add_edge(event_list.index(i[0]['Id']), event_list.index(i[1]['Id']))
    return G

def create_graph_from_events2(events):
    G = nx.Graph()
    event_list = [event['Id'] for event in events]
    G.add_nodes_from([(event['Id'], {'color' : None}) for event in events])

    comb = combinations(events, 2)
    for i in comb:
        res0 = set(i[0]['Resources'])
        res1 = i[1]['Resources']
        intersection = [value for value in res0 if value in res1]
        if intersection:
            G.add_edge(i[0]['Id'], i[1]['Id'])
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
    school = "Den"

    # --------------------------
    # Parse XML file
    # --------------------------
    events = parseXML('dataset/den-smallschool.xml')

    # --------------------------
    #  Preparing Conflict Graph
    # --------------------------
    G = create_graph_from_events(events)
    G_tuple = create_graph_from_events2(events)

    # ---------------------------
    # Verifying Graph consistency
    #----------------------------
    print("----------------------------")
    print("Verifying Graph consistency")
    print("Nodes = ", G.nodes)

    for i in G.nodes:
        print("\nNode", i,"Color", G.nodes[i]['color'])
        color_and_neighbour = [(neighbour, G.nodes[neighbour]['color']) for neighbour in G[i]]
        print("Neighbours | Color")
        for pair in color_and_neighbour:
            print(pair)
    print("----------------------------")
    print("Verifying Graph consistency")
    print("Nodes = ", G_tuple.nodes)

    for i in G_tuple.nodes:
        print("\nNode", i,"Color", G_tuple.nodes[i]['color'])
        color_and_neighbour = [(neighbour, G_tuple.nodes[neighbour]['color']) for neighbour in G_tuple[i]]
        print("Neighbours | Color")
        for pair in color_and_neighbour:
            print(pair)


    num_colors = 5
    print("\nNumber of colors", num_colors)

    # Probabilities list to be used at the end of the program 
    states = []
    color = [None] * G.number_of_nodes()
    # print all k–colorable configurations of the graph
    kColorable(G, color, num_colors, 0, G.number_of_nodes(), states)

    print("Number of states:", len(states))
    #pp.pprint(states)
    
    print("\nRunning Random Distribution")
    
    # -----------------------------------------------------
    # Evaluate the data from limited number of Measurements
    # -----------------------------------------------------
    for it in range(10):
        print("Iteration number:", it)
        measurement_number = 10000
        measurement = np.random.choice(states, measurement_number)
        unique, repet = np.unique(measurement, return_counts=True)
        analysis = dict(zip(unique, repet))

        avr_C       = 0
        hist        = {}

        for sample in list(analysis.keys()):
            coloring = [int(num) for num in sample.split(" ")]
            color_graph_from_coloring(G_tuple, coloring)

            #fx = cost_function_den_25pts(G_tuple)
            fx = cost_function_den_4pts(G_tuple)

            # compute the expectation value and energy distribution
            avr_C     = avr_C    + analysis[sample]*fx
            hist[str(round(fx))] = hist.get(str(round(fx)),0) + analysis[sample]

        
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
        print("----------------------------")

if __name__ == '__main__':
    main()
