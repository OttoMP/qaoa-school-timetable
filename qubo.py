import xhsttparser
import itertools
from pyqubo import OneHotEncInteger, LogEncInteger, NotConst, XorConst
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
            H = constraint.weight*sum(all_events_as_variables)
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
                        #print("Target", target)
                        #print("Target found", instance.events.event[target].name.strip())
                        for variable in all_events_as_variables:
                            if instance.events.event[target].name.strip() == variable.label:
                                #print("Variable Found")
                                for time in constraint.times:
                                    index = available_times.index(time)
                                    #print(f"Index time {bin(index)}") 
                                    exp = NotConst(variable, index, constraint.name)
                                    #print("Constraint", exp)
                                    H += constraint.weight*exp
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
                resources0 = pair[0].resources
                resources1 = pair[1].resources
                conflict = set(resources0+resources1)
                if conflict:
                    pair_names = [event.name.strip() for event in pair]
                    variables = [v for v in all_events_as_variables if v.label in pair_names]
                    #print("Clashing Resource")
                    #print(variables)
                    exp = XorConst(variables[0], variables[1], 1, constraint.name)
                    H += exp

            #print("-------------")

    #print("\n---------------------\n")
    model = H.compile()
    
    return model