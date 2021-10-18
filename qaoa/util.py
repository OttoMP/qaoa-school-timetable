import networkx as nx
from itertools import combinations

def create_graph_from_events(events):
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

def color_graph_from_coloring(graph, coloring):
    for index,node in enumerate(graph.nodes):
        graph.nodes[node]['color'] = coloring[index]
    return

def color_graph_greedy(G):
    pair = None, G.number_of_nodes(), 0
    it = 0
    for i in range (1, 10000):
        color_by_node, colors = color_graph_greedy_aux(G, 0.7)
        if pair[1] > len(colors):
            pair = color_by_node, len(colors), it
        it+= 1
    # Coloring Graph
    for key, value in pair[0].items(): 
        G.nodes[key]['color'] = value

def color_graph_greedy_aux(G, alpha):
    
    n = G.number_of_nodes()
    
    #lista de cores
    colors = []
    
    #vertices ordenadas por degree decrescente
    nodes_ordered = sorted(G.degree, key=lambda x: x[1], reverse=True)
    
    
    #lista de cores proibidas
    forbidden_colors = {}
    for (node, d) in nodes_ordered:
        forbidden_colors[node]=[]
        
    color_by_node = {}
    while len(nodes_ordered) > 0:
        node = None
        i = 0
        while i < len(nodes_ordered):
            if random.random() <=  alpha or i == len(nodes_ordered)-1:
                (node, d) = nodes_ordered.pop(i)
                break
            i+=1              
        
        p_colors = list(set.difference(set(colors), set(forbidden_colors[node]))) 
        c = 0
        if len(p_colors) > 0:
            c = p_colors[0]
            color_by_node[node] = c
        else:
            c = len(colors)
            colors.append(c)
            color_by_node[node] = c
        #proibe cor para adjacentes
        for adj in G.neighbors(node):
            forbidden_colors[adj].append(c)

    return color_by_node, colors

def color_graph_from_num(graph, num_color):
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
