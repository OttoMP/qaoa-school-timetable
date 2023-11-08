import warnings
import xhsttparser
import itertools
from pyqubo import (
    And,
    Xor,
    Not,
    LogEncInteger,
    NotConst,
)


class XhsttToQubo:
    def __init__(self, instance: xhsttparser.Instance) -> None:
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

        self.num_time_slots = len(self.available_times)

        self.events_as_variables = {
            event: LogEncInteger(event, (0, self.num_time_slots - 1))
            for event in instance.events.event.keys()
        }

    def compile(self):
        h = 0
        for constraint in self.instance.constraints:
            if isinstance(constraint, xhsttparser.AssignTimeConstraint):
                h += self.assign_time_constraint(constraint)
            # elif isinstance(constraint, xhsttparser.PreferTimesConstraint):
            #    h += self.prefer_times_constraint(constraint)
            elif isinstance(constraint, xhsttparser.AvoidClashesConstraint):
                h += self.avoid_clash_constraint(constraint)
            else:
                warnings.warn(f"Constraint {constraint} not supported")
        return h.compile()

    def assign_time_constraint(self, constraint: xhsttparser.AssignTimeConstraint):
        h = 0
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
                var_a = self.events_as_variables[event_name_a]
                var_b = self.events_as_variables[event_name_b]

                h += -((var_a - var_b) ** 2)

        return h
