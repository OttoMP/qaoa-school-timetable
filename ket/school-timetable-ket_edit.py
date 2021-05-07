#import math tools
import numpy as np

# We import the tools to handle general Graphs
import networkx as nx

# Import miscellaneous tools
import pprint as pp
import time
import progressbar

# We import plotting tools
import pandas as pd
import matplotlib.pyplot as plt
from   matplotlib import cm
from   matplotlib.ticker import LinearLocator, FormatStrFormatter
from qiskit.visualization import plot_histogram

# Import tools for running QAOA
import random
from scipy.optimize import minimize, fmin, Bounds
from ket import *
from ket.lib import swap, within

# Parallelization tools
import ray
import multiprocessing

#num_cores = multiprocessing.cpu_count()
num_cores = 2
ray.init(num_cpus=num_cores)

def show_figure(fig):
    new_fig = plt.figure()
    new_mngr = new_fig.canvas.manager
    new_mngr.canvas.figure = fig
    fig.set_canvas(new_mngr.canvas)
    plt.show(fig)

# Compute the value of the cost function
# Minimize the number of USED colors
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

# Compute the value of the cost function
def cost_function_C(x,G,num_colors):
    C = 0
    for color in range(num_colors):
        color_used = 0
        for node in G.nodes():
            if x[node*num_colors + color] == 1:
                color_used = 1
                break
        C += color_used

    return C

def save_csv(data, nome_csv):
    data_points = pd.DataFrame(data, columns=['Expected Value', 'p', 'Graph Number'])
    data_points.to_csv(nome_csv, mode='a', header=False)

    return

def color_graph_num(graph, num_color, root):
    color_index = 0
    not_allowed_color = []

    # Mark all the vertices as not visited
    visited = {x: False for x in graph.nodes}

    # Create a queue for BFS
    queue = []

    # Mark the source node as
    # visited and enqueue it
    queue.append(root)
    visited[root] = True

    while queue:
        # Dequeue a vertex from queue and color it
        source = queue.pop(0)
        not_allowed_color = [graph.nodes[neighbour]['color'] for neighbour in graph[source]
                                if (graph.nodes[neighbour]['color'] != None
                                    and graph.nodes[neighbour]['color'] < num_color) ]
        while color_index in not_allowed_color:
            color_index = (color_index+1)%num_color
        graph.nodes[source]['color'] = color_index
        not_allowed_color = []

        # Get all adjacent vertices of the
        # dequeued vertex s. If a adjacent
        # has not been visited, then mark it
        # visited and enqueue it
        for i in graph[source]:
            if visited[i] == False:
                queue.append(i)
                visited[i] = True

    return

def color_graph_coloring(graph, coloring):
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

def create_graph_tuple(nodes, edges):
    G = nx.Graph()
    G.add_nodes_from([(tuple, {'color' : None}) for tuple in nodes])

    for e, row in enumerate(edges):
        for f, column in enumerate(row):
            if column == 1:
               G.add_edge(nodes[e],nodes[f])

    return G

def partial_mixer(qc, neighbour, ancilla, target, beta):
    def outer():
        if neighbour == None:
            x(ancilla)
        else:
            within(lambda : x(neighbour),
                    lambda : ctrl(neighbour, x, ancilla))

    within(lambda : outer(),
            lambda : [within(lambda :[h(target),
                                        ctrl(target[0], x, target[1])],
                            lambda : ctrl(ancilla, rz, 2*beta, target[1])),
                        within(lambda :[rx(-np.pi/2, target),
                                        ctrl(target[0], x, target[1])],
                            lambda : ctrl(ancilla, rz, 2*beta, target[1]))])

def neighbourhood(G, num_colors, node, color):
    neighbour = G[node]
    neighbour_qubit = [color+(num_colors*u) for u in neighbour]

    return neighbour_qubit

# Apply the partial mixer for each pair of colors of each node
def mixer(qc, G, beta, num_nodes, num_colors):
    for u in range(num_nodes):
        for i in range(num_colors):
            neighbours_i = neighbourhood(G, num_colors, u, i)
            for j in range(num_colors):
                if i < j:
                    neighbours_j = neighbourhood(G, num_colors, u, j)
                    neighbours = neighbours_i+neighbours_j

                    if neighbours == []:
                        q_neighbours = None
                    else:
                        q_neighbours = qc[neighbours[0]]
                        for node in neighbours[1:]:
                            q_neighbours = q_neighbours | qc[node]
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

def phase_separator(qc, G, gamma, num_nodes, num_colors):
    qubits2color(qc, num_nodes, num_colors)
    for node in range(num_colors*num_nodes):
        x(qc[node])
    for k in range(num_colors):
        qubits = [k*num_nodes+node for node in range(num_nodes)]
        control = qc[qubits[0]]
        for qub in qubits[1:-1]:
            control = control | qc[qub]
        target = qc[qubits[-1]]
        ctrl(control, rz, 2*gamma, target)
        #cRz(qc, qubits, 2*gamma)
    for node in range(num_colors*num_nodes):
        x(qc[node])
    color2qubits(qc, num_nodes, num_colors)

def qaoa_min_graph_coloring(p, G, num_colors, gamma, beta0, beta):
    num_nodes = G.number_of_nodes()
    qc = quant((num_nodes*num_colors) + num_nodes)

    #print(len(qc))
    # Initial state preparation
    coloring = [G.nodes[node]['color'] for node in G.nodes]
    for i, color in enumerate(coloring):
        x(qc[(i*num_colors)+color])

    # Alternate application of operators
    mixer(qc, G, beta0, num_nodes, num_colors) # Mixer 0
    for step in range(p):
        phase_separator(qc, G, gamma[step], num_nodes, num_colors)
        mixer(qc, G, beta[step], num_nodes, num_colors)

    #print(report())
    # Measurement
    #result = measure(qc).get()
    return dump(qc)

def qaoa(par, p, shots, G, G_tuple, num_colors, students_list):
    # QAOA parameters
    beta0, par= par[0], par[1:]
    middle = int(len(par)/2)
    gamma = par[:middle]
    beta = par[middle:]

    num_nodes = G.number_of_nodes()
    counts = {}

    # run on local simulator                                                 #############
    result = qaoa_min_graph_coloring(p, G, num_colors, gamma, beta0, beta)   #############
    for i in result.get_states():                                            #############
        binary = np.binary_repr(i, width=(num_nodes*num_colors)+num_nodes)   #############
        counts[binary] = int(2**20*result.probability(i))                    #############

    #pp.pprint(counts)
    #print("==Result==")
    # Evaluate the data from the simulator
    avr_C       = 0
    min_C       = [0, G.number_of_nodes()+1]

    for sample in list(counts.keys()):
        if counts[sample] > 0:
            # use sampled bit string x to compute C(x)
            x         = [int(num) for num in list(sample)]
            tmp_eng   = cost_function_timetable(x,G_tuple, num_colors, students_list)

            # compute the expectation value and energy distribution
            avr_C     = avr_C    + counts[sample]*tmp_eng

            # save best bit string
            if( min_C[1] > tmp_eng):
                min_C[0] = sample
                min_C[1] = tmp_eng

    #M1_sampled   = avr_C/shots
    M1_sampled   = avr_C/sum(counts.values())
    #print(M1_sampled)
    return M1_sampled

@ray.remote
def minimization_process(p, shots, G, G_tuple, num_colors, students_list):
    gamma = [random.uniform(0, 2*np.pi) for _ in range(p)]
    beta0 =  [random.uniform(0, np.pi)]
    beta  = [random.uniform(0, np.pi) for _ in range(p)]
    qaoa_par = beta0+gamma+beta
    qaoa_args = p, shots, G, G_tuple, num_colors, students_list
    print("\nMinimizing function\n")
    res = minimize(qaoa, qaoa_par, args=qaoa_args, method='Nelder-Mead',
            options={'maxiter': 300, 'xatol': 0.1, 'fatol': 0.01, 'disp': True, 'adaptive':True})
    print(res)
    return [res['fun'], p, res['x']]

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

    G = create_graph(lectures, lectureConflict)
    G_tuple = create_graph_tuple(lectures, lectureConflict)

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


    print("Running QAOA")
    p = 1
    shots = 100

    # Parallel task using ray
    Mp1_sampled = ray.get([minimization_process.remote(p, shots, G, G_tuple, num_colors, students_list) for iteration in progressbar.progressbar(range(2))])
    print("Mp1_sampled")
    #pp.pprint(Mp1_sampled)

    '''
    #archive_name = "results/p1.csv"
    Mp1_sampled = []
    for iteration in progressbar.progressbar(range(1)):
        gamma = [random.uniform(0, 2*np.pi) for _ in range(p)]
        beta0 =  [random.uniform(0, np.pi)]
        beta  = [random.uniform(0, np.pi) for _ in range(p)]
        qaoa_par = beta0+gamma+beta
        qaoa_args = p, shots, G, G_tuple, num_colors, students_list
        print("\nMinimizing function\n")
        res = minimize(qaoa, qaoa_par, args=qaoa_args, method='Nelder-Mead',
                options={'maxiter': 300, 'xatol': 0.01, 'fatol': 0.01, 'disp': True, 'adaptive':True})
        #print(res)
        Mp1_sampled.append([res['fun'], p, res['x']])
        #save_csv(Mp1_sampled, archive_name)

    print("Mp1_sampled")
    #pp.pprint(Mp1_sampled)
    '''

    final_answer = min(Mp1_sampled, key=lambda x: x[0])
    beta0 = final_answer[2][0]
    gamma = final_answer[2][1:p+1]
    beta  = final_answer[2][p+1:]

    print("Using Following parameters:")
    print("Beta 0:", beta0)
    print("Gamma:", gamma)
    print("Beta:", beta)
    print("\n")

    num_nodes = G.number_of_nodes()
    counts = {}

    # run on local simulator                                                 #############
    result = qaoa_min_graph_coloring(p, G, num_colors, gamma, beta0, beta)   #############
    for i in result.get_states():                                            #############
        binary = np.binary_repr(i, width=(num_nodes*num_colors)+num_nodes)   #############
        prob = int(2**20*result.probability(i))
        if prob > 0:
            counts[binary] = prob                    #############

    #pp.pprint(counts)

    print("==Result==")
    #print("Average Value ", final_answer[0])
        # Evaluate the data from the simulator
    avr_C       = 0
    min_C       = [0, G.number_of_nodes()+1]
    hist        = {}

    for k in range(num_colors+1):
        hist[str(k)] = hist.get(str(k),0)

    for sample in list(counts.keys()):
        if counts[sample] > 0:
            # use sampled bit string x to compute C(x)
            x         = [int(num) for num in list(sample)]
            tmp_eng   = cost_function_timetable(x,G_tuple, num_colors, students_list)

            # compute the expectation value and energy distribution
            avr_C     = avr_C    + counts[sample]*tmp_eng
            hist[str(round(tmp_eng))] = hist.get(str(round(tmp_eng)),0) + counts[sample]

            # save best bit string
            if( min_C[1] > tmp_eng):
                min_C[0] = sample
                min_C[1] = tmp_eng

    #M1_sampled   = avr_C/shots
    total_counts = sum(counts.values())
    print("Total Counts", total_counts)
    M1_sampled   = avr_C/total_counts
    print("Expected Value = ", M1_sampled)

    print('\n --- SIMULATION RESULTS ---\n')
    max_counts = max(counts, key=lambda key: counts[key])

    print("Best result found: ", min_C[0])
    print("Number of times result showed: ", counts[min_C[0]])
    print("Percentage of times result showed: ", (counts[min_C[0]]/total_counts)*100)
    print("Objective function value: ", min_C[1])
    print()
    print("Most commom result found: ", max_counts)
    print("Number of times result showed: ", counts[max_counts])
    print("Percentage of times result showed: ", (counts[max_counts]/total_counts)*100)
    max_value = cost_function_timetable(max_counts, G_tuple, num_colors, students_list)
    print("Objective function value: ", max_value)


    list_qubits = max_counts
    for i in range(len(G)):
        for pos, char in enumerate(list_qubits[i*num_colors:(i*num_colors+num_colors)]):
            if int(char):
                # color = pos
                maximum_coloring.append(pos)

    print("\nMaximum Coloring",maximum_coloring)
    print("\nMaximum Coloring Qudits values")
    for i in range(len(G_tuple)):
        print(list_qubits[i*num_colors:(i*num_colors+num_colors)])

    print("\nNew Graph information")
    print("\nDegree of each node", degree)
    #print("\nNumber of colors", num_colors)
    color_graph_coloring(G_tuple, maximum_coloring)
    for i in G_tuple.nodes:
        print("\nNode",i,"Color", G_tuple.nodes[i]['color'])
        neighbours = [G_tuple.nodes[neighbour]['color'] for neighbour in G_tuple[i]]
        print("Neighbours Colors", neighbours)



    list_qubits = min_C[0]
    for i in range(len(G)):
        for pos, char in enumerate(list_qubits[i*num_colors:(i*num_colors+num_colors)]):
            if int(char):
                # color = pos
                best_coloring.append(pos)

    print("\nBest Coloring",best_coloring)
    print("\nBest Coloring Qudits values")
    for i in range(len(G_tuple)):
        print(list_qubits[i*num_colors:(i*num_colors+num_colors)])

    #print('\nThe approximate solution is x* = %s with C(x*) = %d' % (min_C[0],min_C[1]))
    #print('The number of times this solution showed was: %d \n' %(counts[min_C[0]]))
    #print('The sampled mean value is Mp_sampled = %.02f' % (M1_sampled))

    print("\nNew Graph information")
    print("\nDegree of each node", degree)
    #print("\nNumber of colors", num_colors)
    color_graph_coloring(G_tuple, best_coloring)
    for i in G_tuple.nodes:
        print("\nNode",i,"Color", G_tuple.nodes[i]['color'])
        neighbours = [G_tuple.nodes[neighbour]['color'] for neighbour in G_tuple[i]]
        print("Neighbours Colors", neighbours)

    #coloring = [G.nodes[node]['color'] for node in G.nodes]
    #print("\nFinal coloring", final_coloring)

    '''
    print("Histogram", hist)
    hist_max = sum(counts.values())
    hist_3 = (hist['3']/hist_max)*100
    hist_2 = (hist['2']/hist_max)*100
    hist_0 = (hist['0']/hist_max)*100
    hist_1 = (hist['1']/hist_max)*100
    hist_4 = (hist['4']/hist_max)*100
    hist_5 = (hist['5']/hist_max)*100
    print("Histogram 3", hist_3)
    print("Histogram 2", hist_2)
    print("Histogram 1", hist_1)
    print("Histogram 0", hist_0)
    print("Histogram 4", hist_4)
    print("Histogram 5", hist_5)
    #print('The cost function is distributed as: \n')
    #plot_histogram(hist,figsize = (8,6),bar_labels = False)
    #plt.savefig("histogram.pdf")
    '''

if __name__ == '__main__':
    main()