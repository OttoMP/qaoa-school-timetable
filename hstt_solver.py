import os, sys
from xhsttparser import XHSTT
from qubo import xhstt_to_qubo
from pprint import pprint

def main():
    print("Starting program\n")

    # Root directory for dataset
    # --------------------------
    # School Instances
    # --------------------------
    #path = os.path.abspath(os.path.join(os.getcwd(), 'dataset', sys.argv[1]+".xml"))
    path = sys.argv[1]

    # --------------------------
    # Parse XML file
    # --------------------------
    #events = parseXML('dataset/den-smallschool.xml')
    xhstt = XHSTT(path)

    instance = next(iter(xhstt.instances.values()))
    
    #print(instance)
    #print(instance.times, end="\n\n")
    #print(instance.resources, end="\n\n")
    #print(instance.events, end="\n\n")
    #print(instance.constraints)

    # ------------------------------
    #  Preparing Problem Formulation
    # ------------------------------

    qubo, offset = xhstt_to_qubo(instance)
    pprint(qubo)
    pprint(offset)

    #----------------------------
    # Starting QAOA
    #----------------------------
    print("----------------------------")
    print("Running QAOA")
    #number_of_qubits = num_nodes*num_colors+num_nodes
    number_of_qubits = len(qubo)
    print(f"Necessary number of qubits: {number_of_qubits}")
    # QAOA parameter
    final_depth = 8
    print(f"Final value of p reached: {final_depth}")

    # Minimizing Example DEN
    #minimization_process_cobyla(goal_p, G, num_colors, school, cost_function_den_4pts)

    print("Program End")
    print("----------------------------")

if __name__ == '__main__':
    main()
