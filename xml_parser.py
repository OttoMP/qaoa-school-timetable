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
        lect['Name'] = item.find('Name').text

        resources = []
        for resource in item.findall('Resources/Resource'):
            resources.append(resource.attrib['Reference'])
            lect['Resources'] = resources

        lectures.append(lect)

    return lectures

def create_graph(events):
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

def main():
    # parse xml file
    #events = parseXML('dataset/den-smallschool.xml')
    events = parseXML('dataset/bra-instance01.xml')

    pp.pprint(events)
    G = create_graph(events)

    # Graph Information
    print("\nGraph information")

    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nPre-coloring", coloring)

    degree = [deg for (node, deg) in G.degree()]
    print("\nDegree of each node", degree)

    num_colors = 5
    print("\nNumber of colors", num_colors)

    node_list = list(G.nodes)
    color_graph_num(G, num_colors, node_list[0])
    #color_graph_coloring(G, initial_coloring)

    for i in G.nodes:
        print("\nNode",i,"Color", G.nodes[i]['color'])
        neighbours = [G.nodes[neighbour]['color'] for neighbour in G[i]]
        print("Neighbours Colors", neighbours)

    coloring = [G.nodes[node]['color'] for node in G.nodes]
    print("\nInitial coloring", coloring)

    #nx.draw(G, with_labels=True, font_weight='bold')
    #plt.show()


if __name__ == "__main__":

    # calling main function
    main()