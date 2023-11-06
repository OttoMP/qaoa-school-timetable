import xhsttparser
import itertools
from math import prod
from functools import reduce
from pyqubo import And, Or, Xor, Not, OneHotEncInteger, LogEncInteger, NotConst, XorConst, OrConst, AndConst, Constraint
from pprint import pprint

def xhstt_to_qubo(instance):
    available_times = []
    #print(f"Instance Times {instance.times}")
    for time_group in instance.times.time_groups:
        #print(f"Time_groups {instance.times.time_groups[time_group]}")
        for time in instance.times.time_groups[time_group]:
            available_times.append(time)
            #print(f"Time {time}")
    
    existing_timeslots = len(available_times)
    #print(f"Existing Timeslots {existing_timeslots}")
    #print(f"Available Times {available_times}")
    #print("\n---------------------\n")

    all_events_as_variables = []
    #pprint(f"Instance Events {instance.events}")
    #print("\n---------------------\n")
    #pprint("Instance All Events", instance.events.event_groups["gr_AllEvents"])
    for event in instance.events.event_groups["gr_AllEvents"]:
        #print(f"Event {event.name}")
        all_events_as_variables.append(LogEncInteger(event.name.strip(), (0, existing_timeslots-1)))
        #print("---------------------")
    
    #print("All Events", all_events_as_variables)
    #print("\n---------------------\n")

    H = 0
    for constraint in instance.constraints:
        if type(constraint) == xhsttparser.AssignTimeConstraint:
            #H = constraint.weight*sum(all_events_as_variables)
            #print("Element", constraint.element)
            #print("Id", constraint.id)
            #print("Name", constraint.name)
            #print("Required", constraint.required)
            #print("Weight", constraint.weight)
            #print("Cost Function", constraint.cost_function)
            #print("Applies to", constraint.applies_to)

            for group in constraint.applies_to:
                if group == "event_groups":
                    pass
                    #for target in constraint.applies_to[group]:
                        #print(instance.events.event_groups[target])
                elif group == "events":
                    pass
                elif group == "resource_groups":
                    pass
            #print("-------------")
        
        elif type(constraint) == xhsttparser.PreferTimesConstraint:
            #print("Element", constraint.element)
            #print("Id", constraint.id)
            #print("Name", constraint.name)
            #print("Time", constraint.times)
            #print("Time Group", constraint.time_groups)
            #print("Required", constraint.required)
            #print("Weight", constraint.weight)
            #print("Cost Function", constraint.cost_function)
            #print("Applies to", constraint.applies_to)
            
            for group in constraint.applies_to:
                if group == "event_groups":
                    pass
                elif group == "events":
                    for target in constraint.applies_to[group]:
                        #print(f"{target=}")
                        #print(f"Found {instance.events.event[target].name.strip()}")
                        for variable in all_events_as_variables:
                            if instance.events.event[target].name.strip() == variable.label:
                                for time in constraint.times:
                                    #print(f"{time=}")
                                    index = available_times.index(time)
                                    #print(f"Index {bin(index)}")
                                    exp = (variable - index)**2
                                    #print(f"{exp=}")
                                    #H += 2*constraint.weight*exp
                elif group == "resource_groups":
                    pass
        #print("-------------")
        
        elif type(constraint) == xhsttparser.AvoidClashesConstraint:
            #print("Element", constraint.element)
            #print("Id", constraint.id)
            #print("Name", constraint.name)
            #print("Required", constraint.required)
            #print("Weight", constraint.weight)
            #print("Cost Function", constraint.cost_function)
            #print("Applies to", constraint.applies_to)

            pair_order_list = itertools.combinations(instance.events.event_groups["gr_AllEvents"],2)

            for pair in pair_order_list:
                resources0 = [res.element.attrib["Reference"] for res in pair[0].resources]
                resources1 = [res.element.attrib["Reference"] for res in pair[1].resources]
                conflict = [i for i in resources0 if i in resources1]
                #conflict = set(resources0+resources1)
                #print(f"{pair=}")
                #print("Resources")
                #print(f"{resources0=}")
                #print(f"{resources1=}")
                #print(f"{conflict=}")
                if conflict:
                    pair_names = [event.name.strip() for event in pair]
                    print(f"{pair_names=}\n")
                    variables = [v for v in all_events_as_variables if v.label in pair_names]
                    #print(f"{variables=}\n")
                    #print("Clashing Resource")
                    #print(variables)
                    #partials = sum([XorConst(a, b, 0, f'xor_{a}_{b}') for a, b in zip(variables[0].array, variables[1].array)])
                    #print(f"{partials=}")
                    # XOR = (A OR B) AND (NOT (A AND B))
                    #xors = [And(Or(a,b), Not(And(a,b))) for a, b in zip(variables[0].array, variables[1].array)]
                    xors = [Not(Xor(a,b)) for a, b in zip(variables[0].array, variables[1].array)]
                    clashs = And(xors[0], xors[1])
                    #for gate in xors:
                        #print(f"{gate=}\n")
                    #factor = len(variables[0].array)
                    #print(f"{factor=}")

                    exp = NotConst(clashs, 1, f'clash_{pair_names[0]}_{pair_names[1]}')
                    #exp = AndConst(xors[0], xors[1], 0, f'clash_{pair_names[0]}_{pair_names[1]}')
                    #exp = (2-sum([Xor(a,b) for a, b in zip(variables[0].array, variables[1].array)]))
                    #exp = sum([(1-Xor(a,b))**2 for a, b in zip(variables[0].array, variables[1].array)])-2
                    #exp = reduce(lambda a, b : AndConst(a, b, 1, constraint.name), [(XorConst(a, b, 0, 'xor_partial')) for a, b in zip(variables[0].array, variables[1].array)])
                    #exp  = AndConst(partials[0], partials[1], 1, constraint.name)
                    #exp = Constraint((variables[0] - variables[1])**2, f"Clash_{pair_names[0]}_{pair_names[1]}", condition=lambda x: x!=10.0)
                    #exp = sum([(AndConst(a,b,0, f'xor_{a}_{b}')) for a, b in zip(variables[0].array, variables[1].array)])
                    #exp = AndConst(a, b, 1, f'and {pair_names[0]} {pair_names[1]}')
                    #exp = - weight*((variables[0] - variables[1])**2)
                    print(f"{exp=}\n")
                    H += 10*constraint.weight*exp
        
        #print("-------------")

    #print("\n---------------------\n")
    model = H.compile()
    
    return model