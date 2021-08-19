# Import tools for running QAOA
import random
import cma
from ket import *

#import math tools
import numpy as np

# We import the tools to handle general Graphs
import networkx as nx

# Import miscellaneous tools
from scipy.optimize import minimize
import sys, os, psutil, datetime
import pprint as pp
import pandas as pd

def save_csv(data, nome_csv):
    data_points = pd.DataFrame(data, columns=['Expected Value', 'p', 'Beta0|Gamma|Beta'])
    data_points.to_csv(nome_csv, mode='a', header=False)

    return

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

def qaoa(par, p, initial_G, num_colors, students_list):
    # --------------------------
    # Unpacking QAOA parameters
    # --------------------------
    beta0 = par[0]
    new_par = np.delete(par, 0)
    middle = int(len(par)/2)
    gamma = new_par[:middle]
    beta = new_par[middle:]
    num_nodes = initial_G.number_of_nodes()

    # --------------------------
    # Verifying Parameters
    # --------------------------
    #print("Using Following parameters:\nBeta0:", beta0, "\nGamma:", gamma, "\nBeta:", beta)

    # --------------------------
    # Running QAOA on simulator
    # --------------------------
    G = nx.Graph()
    G.add_nodes_from(initial_G)
    G.add_edges_from(initial_G.edges)
    initial_coloring = [initial_G.nodes[node]['color'] for node in initial_G.nodes]
    color_graph_from_coloring(G, initial_coloring)
    
    result = qaoa_min_graph_coloring(p, initial_G, num_nodes, num_colors, beta0, gamma, beta)

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
    min_function_value = [0, np.inf]
    hist        = {}

    for sample in list(counts.keys()):
        if counts[sample] > 0:
            # use sampled bit string x to compute f(x)
            x       = [int(num) for num in list(sample)]

            coloring = []
            for i in range(len(G)):
                for pos, char in enumerate(x[i*num_colors:(i*num_colors+num_colors)]):
                    if int(char):
                        coloring.append(pos)

            color_graph_from_coloring(G, coloring)

            fx = cost_function_timetable(G, num_colors, students_list)

            # compute the expectation value and energy distribution
            avr_function_value = avr_function_value + counts[sample]*fx
            hist[str(round(fx))] = hist.get(str(round(fx)),0) + counts[sample]

            # save best bit string
            if( min_function_value[1] > fx):
                min_function_value[0] = sample
                min_function_value[1] = fx

    expectation_value = avr_function_value/sum(counts.values())

    #pp.pprint(hist)
    #print("Expectation value", expectation_value)
    #print("---")
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

    print("Using Following parameters:")
    print("Beta0:", beta0)
    print("Gamma:", gamma)
    print("Beta:", beta)

    # --------------------------
    # COBYLA Optimization
    # --------------------------
    print("\nMemory Usage", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    print("Minimizing function using COBYLA\n")
    print("Current Time:-", datetime.datetime.now())
    beta0_bounds = [[0, np.pi]]
    beta_bounds = [[0, np.pi]]*p
    gamma_bounds = [[0, 2*np.pi]]*p
    bounds = beta0_bounds+gamma_bounds+beta_bounds
    #construct the bounds in the form of constraints
    cons = []
    for factor in range(len(bounds)):
        lower, upper = bounds[factor]
        l = {'type': 'ineq',
            'fun': lambda x, lb=lower, i=factor: x[i] - lb}
        u = {'type': 'ineq',
            'fun': lambda x, ub=upper, i=factor: ub - x[i]}
        cons.append(l)
        cons.append(u)
    res = minimize(qaoa, qaoa_par, args=qaoa_args, method='COBYLA',
            constraints=cons, options={'maxiter': 300, 'disp': True})
    print(res)
    print("Current Time:-", datetime.datetime.now())
    print("Memory Usage", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    print("Saving Final Results")
    save_csv([[res['fun'], p, res['x']]], "results/"+school+"p"+str(p)+".csv" )

    # --------------------------
    # Powell Optimization
    # --------------------------
    print("\nMemory Usage", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    print("Minimizing function using Powell\n")
    beta0_bounds = [(0, np.pi)]
    beta_bounds = [(0, np.pi)]*p
    gamma_bounds = [(0, 2*np.pi)]*p
    bounds = beta0_bounds+gamma_bounds+beta_bounds
    print("Current Time:-", datetime.datetime.now())
    res = minimize(qaoa, qaoa_par, args=qaoa_args, method='Powell',
            bounds = bounds, options={'maxiter': 300, 'disp': True})
    print(res)
    print("Current Time:-", datetime.datetime.now())
    print("Memory Usage", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    print("Saving Final Results")
    save_csv([[res['fun'], p, res['x']]], "results/"+school+"p"+str(p)+".csv" )

    # --------------------------
    # CMA-ES Optimization 
    # --------------------------
    print("\nMemory Usage", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    lower_bounds = [0] * ((2*p)+1)
    upper_bounds_beta0 = [np.pi]
    upper_bounds_gamma = [2*np.pi]*p
    upper_bounds_beta  = [np.pi]*p
    upper_bounds = upper_bounds_beta0+upper_bounds_gamma+upper_bounds_beta 
    opts = {'bounds' : [lower_bounds, upper_bounds], 'maxiter': 300, } #'maxfevals': 300}
    sigma0 = 0.3*(2*np.pi)
    print("Initial Step =", sigma0)
    
    es = cma.CMAEvolutionStrategy(qaoa_par, sigma0, opts)
    while not es.stop():
        solutions = es.ask()
        function_values = [qaoa(s, p, G, num_colors, students_list) for s in solutions]
        es.tell(solutions, function_values)
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
    
    print("Memory Usage", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    print("Saving Final Results")
    save_csv([[res[1], p, res[0]]], "results/"+school+"p"+str(p)+".csv" )

def color_graph_from_coloring(graph, coloring):
    for index,node in enumerate(graph.nodes):
        graph.nodes[node]['color'] = coloring[index]

    return

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
    G = create_graph_from_list(lectures, lectureConflict)

    return G, students_list

def main():
    print("Starting program\n")

    # --------------------------
    # School Instances
    # --------------------------
    school = "CEC"

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

    initial_function_value = cost_function_timetable(G, num_colors, students_list)
    print("\nInitial Function Value", initial_function_value)
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
            
    #----------------------------
    # Starting QAOA
    #----------------------------
    print("----------------------------")
    print("Running QAOA")
    num_nodes = G.number_of_nodes()
    number_of_qubits = num_nodes*num_colors+num_nodes
    print("Necessary number of qubits: ", number_of_qubits)

    # QAOA parameter
    p = int(sys.argv[1])

    # Minimizing Example CEC
    print("Running minimization process with p-value", p)
    minimization_process(p, G, num_colors, school, students_list)

    print("Program End")
    print("----------------------------")

if __name__ == '__main__':
    main()