#import math tools
import numpy as np

# We import the tools to handle general Graphs
import networkx as nx

# Import miscellaneous tools
import pprint as pp

# We import plotting tools
import matplotlib.pyplot as plt
from   matplotlib import cm
from   matplotlib.ticker import LinearLocator, FormatStrFormatter

# Compute the value of the cost function
def cost_function_C(x,G):

    E = G.edges()
    if( len(x) != len(G.nodes())):
        return np.nan

    C = 0;
    for index in E:
        e1 = index[0]
        e2 = index[1]

        w      = G[e1][e2]['weight']
        C = C + w*x[e1]*(1-x[e2]) + w*x[e2]*(1-x[e1])

    return C

def qaoa(V, E, beta, gamma):
    # prepare the quantum and classical resisters
    QAOA = quant(len(V))
    #QAOA = QuantumCircuit(len(V), len(V))

    # apply the layer of Hadamard gates to all qubits
    h(QAOA)
    #QAOA.h(range(len(V)))
    #QAOA.barrier()

    # apply the Ising type gates with angle gamma along the edges in E
    for edge in E:
        k = edge[0]
        l = edge[1]
        ctrl(QAOA[k], u1, -2*gamma, QAOA[l])
        u1(gamma, QAOA[k])
        u1(gamma, QAOA[l])
        #QAOA.cp(-2*gamma, k, l)
        #QAOA.p(gamma, k)
        #QAOA.p(gamma, l)

    # then apply the single qubit X rotations with angle beta to all qubits
    #QAOA.barrier()
    rx(2*beta, QAOA)
    #QAOA.rx(2*beta, range(len(V)))

    # Finally measure the result in the computational basis
    result = measure(QAOA).get()
    return result

def main():
    print("Starting program")

    print("Generating Graph")
    # Generating the butterfly graph with 5 nodes
    n     = 5
    V     = np.arange(0,n,1)
    E     =[(0,1,1.0),(0,2,1.0),(1,2,1.0),(3,2,1.0),(3,4,1.0),(4,2,1.0)]

    G     = nx.Graph()
    G.add_nodes_from(V)
    G.add_weighted_edges_from(E)

    # Generate plot of the Graph
    colors       = ['r' for node in G.nodes()]
    default_axes = plt.axes(frameon=True)
    pos          = nx.spring_layout(G)

    #nx.draw_networkx(G, node_color=colors, node_size=600, alpha=1, ax=default_axes, pos=pos)

    print("Function Evaluation")
    # Evaluate the function
    step_size   = 0.1;

    a_gamma         = np.arange(0, np.pi, step_size)
    a_beta          = np.arange(0, np.pi, step_size)
    a_gamma, a_beta = np.meshgrid(a_gamma,a_beta)

    F1 = 3-(np.sin(2*a_beta)**2*np.sin(2*a_gamma)**2-0.5*np.sin(4*a_beta)*np.sin(4*a_gamma))*(1+np.cos(4*a_gamma)**2)

    # Grid search for the minimizing variables
    result = np.where(F1 == np.amax(F1))
    a      = list(zip(result[0],result[1]))[0]

    gamma  = a[0]*step_size;
    beta   = a[1]*step_size;

    # Plot the expetation value F1
    fig = plt.figure()
    ax  = fig.gca(projection='3d')

    surf = ax.plot_surface(a_gamma, a_beta, F1, cmap=cm.coolwarm, linewidth=0, antialiased=True)

    ax.set_zlim(1,4)
    ax.zaxis.set_major_locator(LinearLocator(3))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    #plt.show()

    print("Running QAOA")
    shots = 10000
    counts = {}

    # run on local simulator
    for _ in range(shots):
        result = qaoa(V, E, beta, gamma)
        binary = np.binary_repr(result, width=n)
        if binary in counts:
            counts[binary] += 1
        else:
            counts[binary] = 1

    pp.pprint(counts)

    avr_C       = 0
    max_C       = [0,0]

    for sample in list(counts.keys()):

        # use sampled bit string x to compute C(x)
        x         = [int(num) for num in list(sample)]
        tmp_eng   = cost_function_C(x,G)

        # compute the expectation value and energy distribution
        avr_C     = avr_C    + counts[sample]*tmp_eng

        # save best bit string
        if( max_C[1] < tmp_eng):
            max_C[0] = sample
            max_C[1] = tmp_eng

    M1_sampled   = avr_C/shots

    print('\n --- SIMULATION RESULTS ---\n')
    print('The sampled mean value is M1_sampled = %.02f while the true value is M1 = %.02f \n' % (M1_sampled,np.amax(F1)))
    print('The approximate solution is x* = %s with C(x*) = %d \n' % (max_C[0],max_C[1]))
    print('The cost function is distributed as: \n')

if __name__ == '__main__':
    main()