from qaoa import *

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
    # Greedy coloring to be used in cases where a trivial coloring cannot be
    # found
    # -----------------------------------------------------------------
    #color_graph_greedy(G)
    
    # If a suitable coloring can be found without the greedy method use
    # the color_graph_num method
    # -----------------------------------------------------------------
    #num_colors = 5 # Denmark colors
    #color_graph_from_num(G, num_colors)

    # If a initial state was chosen in advance use color_graph_from_coloring
    # ----------------------------------------------------------------------
    # Coloring 23 points 
    coloring =  [1, 0, 2, 3, 1, 2, 1, 2, 3, 0, 0, 2, 0, 3, 1, 3, 0, 1, 0, 3, 2, 2, 1, 2, 3]
    # Optimal Coloring
    #coloring =  [0, 2, 3, 1, 2, 3, 3, 2, 0, 1, 0, 3, 2, 1, 0, 2, 3, 0, 2, 1, 3, 3, 0, 3, 1]
    color_graph_from_coloring(G, coloring)
    
    #coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nInitial coloring", coloring)

    #num_colors = len(set(coloring))
    num_colors = 5
    print("\nNumber of colors", num_colors)
    
    initial_function_value = cost_function_den_4pts(G)
    print("\nInitial Function Value Max 4", initial_function_value)

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
    goal_p = 8

    # Minimizing Example DEN
    minimization_process_cobyla(goal_p, G, num_colors, school, cost_function_den_4pts)

    print("Program End")
    print("----------------------------")

if __name__ == '__main__':
    main()
