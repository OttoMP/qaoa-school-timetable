import xml.dom.minidom
import xml.etree.ElementTree as ET
import pprint as pp
from itertools import combinations
import matplotlib.pyplot as plt

# We import the tools to handle general Graphs
import networkx as nx

def parseXML(xmlfile):
    # create element tree object
    tree = ET.parse(xmlfile)
    # get root element
    root = tree.getroot()

    # create empty list for news items
    lectures = []
    # iterate news items
    for item in root.findall('./Instances/Instance/Events/Event'):
        lect = {}
        lect['Id'] = item.attrib['Id']
        lect['Duration'] = int(item.find('Duration').text)
        #lect['Name'] = item.find('Name').text

        resources = []
        for resource in item.findall('Resources/Resource'):
            resources.append(resource.attrib['Reference'])
            lect['Resources'] = resources

        lectures.append(lect)

    return lectures

def create_graph(events):
    G = nx.Graph()
    
    week = [
            "Mo_1",
            "Mo_2",
            "Mo_3",
            "Tu_1",
            "Tu_2",
            "Tu_3",
            "We_1",
            "We_2",
            "We_3",
            "Th_1",
            "Th_2",
            "Th_3",
            "Fr_1",
            "Fr_2",
            "Fr_3",
    ]
    comb = combinations(week, 2)
    G.add_nodes_from(week)
    G.add_nodes_from([(day, {'color' : c}) for c, day in enumerate(week)])
    for i in comb:
        G.add_edge(i[0], i[1])
    
    # Adding events to graph
    # No clash constraints Hard 1
    G.add_nodes_from([(event['Id'], {'color' : None}) for event in events])
    comb = combinations(events, 2)
    for i in comb:
        res0 = set(i[0]['Resources'])
        res1 = i[1]['Resources']
        intersection = [value for value in res0 if value in res1]
        if intersection:
            G.add_edge(i[0]['Id'], i[1]['Id'])
            
    # Forbidden time constraint Hard 1
    #  T1 We
    G.add_edge("T1-S1", "We_1")
    G.add_edge("T1-S1", "We_2")
    G.add_edge("T1-S1", "We_3")
    G.add_edge("T1-S2", "We_1")
    G.add_edge("T1-S2", "We_2")
    G.add_edge("T1-S2", "We_3")
    G.add_edge("T1-S3", "We_1")
    G.add_edge("T1-S3", "We_2")
    G.add_edge("T1-S3", "We_3")
    #  T2 Fr
    G.add_edge("T2-S1", "Fr_1")
    G.add_edge("T2-S1", "Fr_2")
    G.add_edge("T2-S1", "Fr_3")
    G.add_edge("T2-S2", "Fr_1")
    G.add_edge("T2-S2", "Fr_2")
    G.add_edge("T2-S2", "Fr_3")
    #  T3 Th
    G.add_edge("T3-S1", "Th_1")
    G.add_edge("T3-S1", "Th_2")
    G.add_edge("T3-S1", "Th_3")
    G.add_edge("T3-S2", "Th_1")
    G.add_edge("T3-S2", "Th_2")
    G.add_edge("T3-S2", "Th_3")
    G.add_edge("T3-S3", "Th_1")
    G.add_edge("T3-S3", "Th_2")
    G.add_edge("T3-S3", "Th_3")
    #  T4 Tu
    G.add_edge("T4-S1", "Tu_1")
    G.add_edge("T4-S1", "Tu_2")
    G.add_edge("T4-S1", "Tu_3")
    G.add_edge("T4-S2", "Tu_1")
    G.add_edge("T4-S2", "Tu_2")
    G.add_edge("T4-S2", "Tu_3")
    G.add_edge("T4-S3", "Tu_1")
    G.add_edge("T4-S3", "Tu_2")
    G.add_edge("T4-S3", "Tu_3")
    #  T5 Mo
    G.add_edge("T5-S2", "Mo_1")
    G.add_edge("T5-S2", "Mo_2")
    G.add_edge("T5-S2", "Mo_3")
    G.add_edge("T5-S3", "Mo_1")
    G.add_edge("T5-S3", "Mo_2")
    G.add_edge("T5-S3", "Mo_3")
    #  T6 We
    G.add_edge("T6-S1", "We_1")
    G.add_edge("T6-S1", "We_2")
    G.add_edge("T6-S1", "We_3")
    G.add_edge("T6-S2", "We_1")
    G.add_edge("T6-S2", "We_2")
    G.add_edge("T6-S2", "We_3")
    G.add_edge("T6-S3", "We_1")
    G.add_edge("T6-S3", "We_2")
    G.add_edge("T6-S3", "We_3")
    #  T7 Tu
    G.add_edge("T7-S1", "Tu_1")
    G.add_edge("T7-S1", "Tu_2")
    G.add_edge("T7-S1", "Tu_3")
    G.add_edge("T7-S3", "Tu_1")
    G.add_edge("T7-S3", "Tu_2")
    G.add_edge("T7-S3", "Tu_3")
    #  T8 Th
    G.add_edge("T8-S1", "Th_1")
    G.add_edge("T8-S1", "Th_2")
    G.add_edge("T8-S1", "Th_3")
    G.add_edge("T8-S2", "Th_1")
    G.add_edge("T8-S2", "Th_2")
    G.add_edge("T8-S2", "Th_3")
    G.add_edge("T8-S3", "Th_1")
    G.add_edge("T8-S3", "Th_2")
    G.add_edge("T8-S3", "Th_3")

    return G

def color_graph_num(graph, num_color, root):
    color_index = 0
    not_allowed_color = []

    # Mark all the vertices as not visited
    visited = {x: False for x in graph.nodes}
    for node in graph.nodes:
        if graph.nodes[node]['color'] != None:
            visited[node] = True

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
                                if (graph.nodes[neighbour]['color'] != None) ]
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

def color_graph_from_num(graph, num_color):
    color_index = 0
    node_list = list(graph.nodes)
    not_allowed_color = []

    # Mark all the vertices as not visited
    visited = {x: False for x in graph.nodes}
    for node in graph.nodes:
        if graph.nodes[node]['color'] != None:
            visited[node] = True

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

def main():
    # parse xml file
    #events = parseXML('dataset/den-smallschool.xml')
    instance = parseXML('dataset/bra-instance01.xml')

    pp.pprint(instance)
    print("---------------------------")

    # Constraints not used: 
    #   Prefer Times for duration 2 Hard 1
    #   Spread events max 1 per day constraint Hard 1
    #   Split Events Duration 1 and 2 Hard 1
    #   At least 1 double lesson constraint Soft 1
    #   At least 2 double lessons constraint Soft 1
    '''
    events = []
    for e in instance:
        if e['Duration'] == 3:
            event1 = {'Id': e['Id']+'pt1', 'Resources': e['Resources'], 'Duration': 1}
            event2 = {'Id': e['Id']+'pt2', 'Resources': e['Resources'], 'Duration': 2}
            events.append(event1)
            events.append(event2)
        elif e['Duration'] == 4:
            event1 = {'Id': e['Id']+'pt1', 'Resources': e['Resources'], 'Duration': 2}
            event2 = {'Id': e['Id']+'pt2', 'Resources': e['Resources'], 'Duration': 2}
            events.append(event1)
            events.append(event2)
        elif e['Duration'] == 5:
            event1 = {'Id': e['Id']+'pt1', 'Resources': e['Resources'], 'Duration': 1}
            event2 = {'Id': e['Id']+'pt2', 'Resources': e['Resources'], 'Duration': 2}
            event3 = {'Id': e['Id']+'pt3', 'Resources': e['Resources'], 'Duration': 2}
            events.append(event1)
            events.append(event2)
            events.append(event3)
        else:
            events.append(e)

    print(len(events)) 
    pp.pprint(events)
    '''
    G = create_graph(instance)

    # Graph Information
    print("\nGraph information")

    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nPre-coloring", coloring)

    degree = [deg for (node, deg) in G.degree()]
    print("\nDegree of each node", degree)

    num_colors = 15
    print("\nNumber of colors", num_colors)

    # Assign Time Constraint Hard 1
    node_list = list(G.nodes)
    #color_graph_num(G, num_colors, node_list[15])
    color_graph_from_num(G, num_colors)
    #color_graph_coloring(G, initial_coloring)

    for i in G.nodes:
        print("\nNode",i,"Color", G.nodes[i]['color'])
        neighbours = [G.nodes[neighbour]['color'] for neighbour in G[i]]
        print("Neighbours Colors", neighbours)

    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nInitial coloring", coloring)

    num_nodes = G.number_of_nodes()
    number_of_qubits = num_nodes*num_colors+num_nodes
    print("Necessary number of qubits: ", number_of_qubits)
    
    #nx.draw(G, with_labels=True, font_weight='bold')
    #plt.show()

if __name__ == "__main__":

    # calling main function
    main()