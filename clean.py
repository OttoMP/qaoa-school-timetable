# Import tools for running QAOA
import random
import cma
from ket import *

#import math tools
import numpy as np

# We import the tools to handle general Graphs
import networkx as nx

# Import miscellaneous tools
import pprint as pp

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

def qaoa(par, p, G, num_colors, students_list):
    # --------------------------
    # Unpacking QAOA parameters
    # --------------------------
    beta0 = par[0]
    new_par = np.delete(par, 0)
    middle = int(len(par)/2)
    gamma = new_par[:middle]
    beta = new_par[middle:]
    num_nodes = G.number_of_nodes()

    # --------------------------
    # Verifying Parameters
    # --------------------------
    print("Using Following parameters:\nBeta0:", beta0, "\nGamma:", gamma, "\nBeta:", beta)

    # --------------------------
    # Running QAOA on simulator
    # --------------------------
    # run on local simulator
    result = qaoa_min_graph_coloring(p, G, num_nodes, num_colors, beta0, gamma, beta)

    #print("Number of States", len(result.get_states()))
    #print("State Vector", result.show('b6:b6:b6:b6:b6:b6'))

    # --------------------------
    # Counting resulting states
    # --------------------------
    counts = {} # Dictionary for keeping the results of the simulation
    for i in result.get_states():
        binary = np.binary_repr(i, width=(num_nodes*num_colors)+num_nodes)
        counts[binary] = int(2**20*result.probability(i))

    # --------------------------
    # Evaluate the data from the simulator
    # --------------------------
    avr_function_value = 0
    min_function_value = [0, np.inf]

    for sample in list(counts.keys()):
        if counts[sample_state] > 0:
            # use sampled bit string x to compute f(x)
            x       = [int(num) for num in list(sample_state)]
            fx = cost_function_timetable(x, G, num_colors, students_list)

            # compute the expectation value and energy distribution
            avr_function_value = avr_function_value + counts[sample_state]*fx

            # save best bit string
            if( min_function_value[1] > fx):
                min_function_value[0] = sample_state
                min_function_value[1] = fx

    expectation_value = avr_function_value/sum(counts.values())

    return expectation_value

def minimization_process(p, G, num_colors, school, students_list):
    # --------------------------
    # Initializing QAOA Parameters 
    # --------------------------
    gamma = [random.uniform(0, 2*np.pi) for _ in range(p)]
    beta0 = random.uniform(0, np.pi)
    beta  = [random.uniform(0, np.pi) for _ in range(p)]
    qaoa_par = [beta0]+gamma+beta
    qaoa_args = p, G, num_colors, students_list

    '''
    # --------------------------
    # Nelder-Mead Optimization
    # --------------------------
    data = []
    res = minimize(qaoa, qaoa_par, args=qaoa_args, method='Nelder-Mead',
            options={'maxiter': 1, 'disp': True, 'adaptive':True})
    print(res)

    data.append([res['fun'], p, res['x']])
    save_csv(data, "results/"+school+"p"+str(p)+".csv" )
    '''
    # --------------------------
    # CMA-ES Optimization 
    # --------------------------
    lower_bounds = [0] * ((2*p)+1)
    upper_bounds_beta0 = [np.pi]
    upper_bounds_gamma = [2*np.pi]*p
    upper_bounds_beta  = [np.pi]*p
    upper_bounds = upper_bounds_beta0+upper_bounds_gamma+upper_bounds_beta 
    opts = {'bounds' : [lower_bounds, upper_bounds], 'maxiter': 2, } #'maxfevals': 300}
    sigma0 = 2
    
    es = cma.CMAEvolutionStrategy(qaoa_par, sigma0, opts)
    while not es.stop():
        solutions = es.ask()
        print("Solutions", solutions)
        es.tell(solutions, [qaoa(s, p, G, num_colors, students_list) for s in solutions])
        res = es.result
        #print("Saving Results")
        save_csv([[res[1], p, res[0]]], "results/"+school+"p"+str(p)+".csv" )
        es.disp()
    print("---------------------------")
    es.result_pretty()
    res = es.result
    print("---------------------------")
    print("Optimal Result", res[0])
    print("Respective Function Value", res[1])
    print("Respective Function Evaluations", res[2])
    print("Overall Function Evaluations", res[3])
    print("Overall Iterations", res[4])
    print("Mean Result", res[5])
    print("Standard Deviation Final Sample", res[6])
    
    print("Saving Final Results")
    save_csv([[res[1], p, res[0]]], "results/"+school+"p"+str(p)+".csv" )

def color_graph_from_coloring(graph, coloring):
    for index,node in enumerate(graph.nodes):
        graph.nodes[node]['color'] = coloring[index]

    return

def create_graphv_from_list(nodes, edges):
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
    G = create_graphv_from_list(lectures, lectureConflict)

    return G, students_list

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
    G, students_list = first_example()
    
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
    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nInitial coloring", coloring)

    # ---------------------------
    # Verifying Graph consistency
    #----------------------------
    print("----------------------------")
    print("Verifying Graph consistency")
    for i in G.nodes:
        print("\nNode",i,"Color", G.nodes[i]['color'])
        color_and_neighbour = [(neighbour, G.nodes[neighbour]['color']) for neighbour in G[i]]
        print("Neighbours | Color")
        for pair in color_and_neighbour:
            print(pair)
            
    #----------------------------
    # Starting QAOA
    #----------------------------
    print("----------------------------")
    print("Running QAOA")
    num_nodes = G.number_of_nodes()
    number_of_qubits = num_nodes*num_colors+num_nodes
    print("Necessary number of qubits: ", number_of_qubits)

    # QAOA parameter
    p = 1

    # Minimizing Example CEC
    print("Running minimization process")
    minimization_process(p, G, num_colors, school, students_list)

if __name__ == '__main__':
    main()