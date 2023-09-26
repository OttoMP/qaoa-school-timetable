import os, sys
import neal
#from dwave.samplers import DWaveSampler
from xhsttparser import XHSTT
from qubo import xhstt_to_qubo
from pprint import pprint
from qat.core import Observable, Term

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

    model = xhstt_to_qubo(instance)
    qubo, offset = model.to_qubo()
    #pprint(qubo)
    #pprint(offset)

    #----------------------------
    # Starting QAOA
    #----------------------------
    print("----------------------------")
    print("Running QAOA")
    number_of_qubits = len(model.variables)
    print(f"Necessary number of qubits: {number_of_qubits}")
    # QAOA parameter
    final_depth = 8
    #print(f"Final value of p reached: {final_depth}")

    single_values = {}
    multiple_values = {}
    indexes = model.variables[::-1]
    for k,v in qubo.items():
        if k[0] == k[1]:
            single_values[indexes.index(k[0])] = v
        else:
            multiple_values[indexes.index(k[0]), indexes.index(k[1])] = v

    hamiltonian_qubo = Observable(number_of_qubits,
                            pauli_terms=
                            [Term(single_values[x], "Z", [x]) for x in single_values]+
                            [Term(multiple_values[x], "ZZ", [x[0],x[1]]) for x in multiple_values],
                            constant_coeff=offset
                            )

    print(hamiltonian_qubo)


    sampler = neal.SimulatedAnnealingSampler()
    bqm = model.to_bqm()
    sampleset = sampler.sample(bqm, num_reads=10)
    decoded_samples = model.decode_sampleset(sampleset)
    best_sample = min(decoded_samples, key=lambda x: x.energy)
    pprint(best_sample.sample)
    
    '''
    bqm = model.to_bqm()
    sa = ExactSolver()
    sampleset = sa.sample(bqm)
    sampleset.to_pandas_dataframe().sort_values("energy")
    '''
    # Minimizing Example DEN
    #minimization_process_cobyla(goal_p, G, num_colors, school, cost_function_den_4pts)

    print("Program End")
    print("----------------------------")

if __name__ == '__main__':
    main()
