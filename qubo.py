import xhsttparser
from pyqubo import LogEncInteger
from pprint import pprint

def xhstt_to_qubo(instance):
    existing_timeslots = 0
    #print(f"Instance Times {instance.times}")
    for time_group in instance.times.time_groups:
        #print(f"Time_groups {instance.times.time_groups[time_group]}")
        for times in instance.times.time_groups[time_group]:
            existing_timeslots += 1
            #print(f"Time {times}")
    
    print(f"Existing Timeslots {existing_timeslots}")
    print("\n---------------------\n")

    all_events_as_variables = []
    #pprint(f"Instance Events {instance.events}")
    #print("\n---------------------\n")
    #pprint("Instance All Events", instance.events.event_groups["gr_AllEvents"])
    for event in instance.events.event_groups["gr_AllEvents"]:
        #print(f"Event {event.name}\n")
        all_events_as_variables.append(LogEncInteger(event.name.strip(), (0, existing_timeslots-1)))
        #print("---------------------")
    
    print("All Events", all_events_as_variables)
    print("\n---------------------\n")

    H = 0
    for constraint in instance.constraints:
        if type(constraint) == xhsttparser.AssignTimeConstraint:
            print("Element", constraint.element)
            print("Id", constraint.id)
            print("Name", constraint.name)
            print("Required", constraint.required)
            print("Weight", constraint.weight)
            print("Cost Function", constraint.cost_function)
            print("Applies to", constraint.applies_to)

            for group in constraint.applies_to:
                print(str(group))
                if group == "event_groups":
                    for target in constraint.applies_to[group]:
                        print(instance.events.event_groups[target])
                    H = sum(all_events_as_variables)
                elif group == "events":
                    pass
                elif group == "resource_groups":
                    pass

                    
            print("-------------")
        elif type(constraint) == xhsttparser.PreferTimesConstraint:
            print("Element", constraint.element)
            print("Id", constraint.id)
            print("Name", constraint.name)
            print("Required", constraint.required)
            print("Weight", constraint.weight)
            print("Cost Function", constraint.cost_function)
            print("Applies to", constraint.applies_to)
            print("-------------")
        elif type(constraint) == xhsttparser.AvoidClashesConstraint:
            print("Element", constraint.element)
            print("Id", constraint.id)
            print("Name", constraint.name)
            print("Required", constraint.required)
            print("Weight", constraint.weight)
            print("Cost Function", constraint.cost_function)
            print("Applies to", constraint.applies_to)
            print("-------------")

    print("\n---------------------\n")
    model = H.compile()

    qubo, offset = model.to_qubo()
    return qubo, offset