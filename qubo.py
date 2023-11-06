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

        self.available_times = []
        for time_group in instance.times.time_groups:
            for time in instance.times.time_groups[time_group]:
                self.available_times.append(time)

        self.existing_time_slots = len(self.available_times)

        self.all_events_as_variables = []
        for event in instance.events.event_groups[
            "gr_AllEvents"
        ]:  # TODO Verify if this cover all instances
            self.all_events_as_variables.append(
                LogEncInteger(event.name.strip(), (0, self.existing_time_slots - 1))
            )

    def compile(self):
        h = 0
        for constraint in self.instance.constraints:
            if isinstance(constraint, xhsttparser.AssignTimeConstraint):
                h += self.assign_time_constraint(constraint)
            # elif isinstance(constraint, xhsttparser.PreferTimesConstraint):
            #     h += self.prefer_times_constraint(constraint)
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
            for variable in self.all_events_as_variables:
                if self.instance.events.event[target].name.strip() == variable.label:
                    for time in constraint.times:
                        index = self.available_times.index(time)
                        exp = (variable - index) ** 2
                        h += 2 * constraint.weight * exp

        return h

    def avoid_clash_constraint(self, constraint: xhsttparser.AvoidClashesConstraint):
        h = 0

        pair_order_list = itertools.combinations(
            self.instance.events.event_groups["gr_AllEvents"], 2
        )

        for pair in pair_order_list:
            resources0 = [res.element.attrib["Reference"] for res in pair[0].resources]
            resources1 = [res.element.attrib["Reference"] for res in pair[1].resources]
            conflict = [i for i in resources0 if i in resources1]

            if conflict:
                pair_names = [event.name.strip() for event in pair]
                # print(f"{pair_names=}\n")
                variables = [
                    v for v in self.all_events_as_variables if v.label in pair_names
                ]

                xors = [
                    Not(Xor(a, b))
                    for a, b in zip(variables[0].array, variables[1].array)
                ]
                clashs = And(xors[0], xors[1])

                exp = NotConst(clashs, 1, f"clash_{pair_names[0]}_{pair_names[1]}")

                h += 10 * constraint.weight * exp
        return h
