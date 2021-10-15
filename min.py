from qaoa import *

def cost_function_min(G):
    C = 0
        
    if G.nodes["Event1"]['color'] != 1:
        C += 1
    if G.nodes["Event2"]['color'] != 2:
        C += 1
    if G.nodes["Event3"]['color'] != 3:
        C += 1
    if G.nodes["Event4"]['color'] != 4:
        C += 1
    if G.nodes["Event5"]['color'] != 0:
        C += 1
    
    return C

def minimal_example():
    import networkx as nx

    nodes = [('Event1', {'color': None}),
             ('Event2', {'color': None}),
             ('Event3', {'color': None}),
             ('Event4', {'color': None}),
             ('Event5', {'color': None}),
             ('Event6', {'color': None}),
             ('Event7', {'color': None}),
    ]
    edges = [('Event1', 'Event2'),
             ('Event1', 'Event3'),
             ('Event1', 'Event4'),
             ('Event1', 'Event5'),
             ('Event2', 'Event3'),
             ('Event2', 'Event4'),
             ('Event2', 'Event5'),
             ('Event3', 'Event4'),
             ('Event3', 'Event5'),
             ('Event4', 'Event5'),
             ('Event1', 'Event6'),
             ('Event2', 'Event6'),
             ('Event3', 'Event6'),
             ('Event4', 'Event6'),
             ('Event5', 'Event6'),
             ('Event1', 'Event7'),
             ('Event2', 'Event7'),
             ('Event3', 'Event7'),
             ('Event4', 'Event7'),
             ('Event5', 'Event7')
    ]

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    return G

def main():
    print("Starting program\n")

    # --------------------------
    # School Instances
    # --------------------------
    school = "Min"

    # --------------------------
    #  Preparing Conflict Graph
    # --------------------------
    G = minimal_example()
    
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
    # Minimal example Coloring
    color_graph_from_coloring(G, [0,1,2,3,4,5,6])
    #color_graph_from_coloring(G, [6,1,2,3,0,5,5])
    
    coloring = [G.nodes[node]['color'] for node in G.nodes]
    num_colors = len(set(coloring))
    print("\nNumber of colors", num_colors)
    print("\nInitial coloring", coloring)

    initial_function_value = cost_function_min(G)
    print("\nInitial Function Value Max 5:", initial_function_value)

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
    goal_p = 8

    # Minimizing Example DEN
    minimization_process_cobyla(goal_p, G, num_colors, school, cost_function_min)
    #minimization_process_cma(goal_p, G, num_colors, school)
    
    print("Program End")
    print("----------------------------")

if __name__ == '__main__':
    main()