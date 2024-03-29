import warnings
import xhsttparser
import itertools
from pyqubo import (
    And,
    Xor,
    Not,
    LogEncInteger,
    OneHotEncInteger,
    NotConst,
    AndConst,
)


class XhsttToQubo:
    def __init__(self, instance: xhsttparser.Instance, OneHotEnc = True) -> None:
        self.instance = instance

        self.available_times = {
            time: index
            for index, time in enumerate(
                [
                    time
                    for time_group in instance.times.time_groups
                    for time in instance.times.time_groups[time_group]
                ]
            )
        }

        self.onehot= OneHotEnc
        self.num_time_slots = len(self.available_times)

        if OneHotEnc:
            self.events_as_variables = {
                event: OneHotEncInteger(event, (0, self.num_time_slots - 1), 5)
                for event in instance.events.event.keys()
            }
        else:
            self.events_as_variables = {
                event: LogEncInteger(event, (0, self.num_time_slots - 1))
                for event in instance.events.event.keys()
            }

    def compile(self):
        h = 0
        for constraint in self.instance.constraints:
            if isinstance(constraint, xhsttparser.AssignTimeConstraint):
                h += self.assign_time_constraint(constraint)
            elif isinstance(constraint, xhsttparser.PreferTimesConstraint):
                h += self.prefer_times_constraint(constraint)
            elif isinstance(constraint, xhsttparser.AvoidClashesConstraint):
                h += self.avoid_clash_constraint(constraint)
            else:
                warnings.warn(f"Constraint {constraint} not supported")
        return h.compile()

    def assign_time_constraint(self, constraint: xhsttparser.AssignTimeConstraint):
        h = 0
        
        if self.onehot:
            for key, value in self.events_as_variables.items():
                h += 1 - sum(value.array)
                for a,b in itertools.combinations(value.array, 2):
                    h += 2*a*b
        else:
            h = sum(self.events_as_variables.values())

        return h

    def prefer_times_constraint(self, constraint: xhsttparser.PreferTimesConstraint):
        h = 0

        for target in constraint.applies_to["events"]:
            for time in constraint.times:
                h += (
                    2
                    * constraint.weight
                    * (self.events_as_variables[target] - self.available_times[time])
                    ** 2
                )

        return h

    def avoid_clash_constraint(self, constraint: xhsttparser.AvoidClashesConstraint):
        h = 0

        event_pairs = itertools.combinations(self.instance.events.event.items(), 2)

        for (event_name_a, event_a), (event_name_b, event_b) in event_pairs:
            resources_a = [res.reference for res in event_a.resources]
            resources_b = [res.reference for res in event_b.resources]
            conflict = [i for i in resources_a if i in resources_b]

            if len(conflict) != 0:
                print(f"{event_name_a=} {event_name_b=}")
                var_a = self.events_as_variables[event_name_a]
                var_b = self.events_as_variables[event_name_b]
                
                if self.onehot:
                    for bit_a, bit_b in zip(var_a.array, var_b.array):
                        h += 5*constraint.weight*(bit_a*bit_b)
                else:
                    xors = [Not(Xor(a,b)) for a, b in zip(var_a.array, var_b.array)]
                    exp = AndConst(xors[0], xors[1], 0, f'clash_{event_name_a}_{event_name_b}')
                    h += 5*constraint.weight*exp

        return h
